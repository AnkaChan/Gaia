#pragma once
#include "aabb.cuh"
#include "GPU_LBVH.cuh"
#include "query.cuh"

#include <CuMatrix/MatrixOps/VectorTypes.h>

#define SQR(x) ((x)*(x))
#define POINT_OVERLAP_EPSILON 1e-6

namespace lbvh
{

    template<typename Real, int32_t PositionBufferStride = 3>
    struct PointAABBGetter 
    {
        __device__ __host__ __forceinline__
            lbvh::aabb<float> operator()(const Real* posBuffer, const IndexType* primitiveBuffer, const IndexType primId) const noexcept
        {
            lbvh::aabb<float> retval;
            //printf("*(primitiveBuffer + primId): %d\n", *(primitiveBuffer + primId));
            const Real* pointBuffer = posBuffer + PositionBufferStride * (*(primitiveBuffer + primId));
            //printf("pointBuffer: %f %f %f\n", pointBuffer[0], pointBuffer[1], pointBuffer[2]);
            retval.upper.x = pointBuffer[0];
            retval.upper.y = pointBuffer[1];
            retval.upper.z = pointBuffer[2];
            retval.lower = retval.upper;
            return retval;
        }
    };


    template<typename Real>
    struct PointCloudClosestPtQueryResult : public QueryDataBase<Real>
    {
        IndexType closestPtGeomId;
        IndexType closestPtPrimId;
        Real distance;
    };

    template<typename Real>
    struct PointOverlapCallBackFunctor {
        __device__ __forceinline__ bool operator()(IndexType obj_idx, IndexType geomId, IndexType primId,
            const detail::basic_device_bvh<Real, true>& bvh, QueryDataBase<Real>* userData) const
        {
            PointCloudClosestPtQueryResult<Real>* queryResult = (PointCloudClosestPtQueryResult<Real>*)userData;

            Real* positionsBuffer = bvh.geometryDataDevice[geomId].positionsBuffer;
            IndexType* primitivesBuffer = bvh.geometryDataDevice[geomId].primitivesBuffer;
            Real* position = positionsBuffer + 3 * *(primitivesBuffer + primId);
            float4 pos = make_float4(position[0], position[1], position[2], 0.f);

            const auto& queryPosition = userData->queryPosition;
            Real dst2 = disVec3Square(queryPosition, pos);
            printf("queryPosition: %f %f %f, pos: %f %f %f, dst2: %f\n", queryPosition.x, queryPosition.y, queryPosition.z, pos.x, pos.y, pos.z, dst2);
            printf("obj_idx: %d, geomId: %d, primId: %d\n", obj_idx, geomId, primId, dst2);
            if (dst2 < 1e-10)
            {
                queryResult->closestPtGeomId = geomId;
                queryResult->closestPtPrimId = primId;
                queryResult->distance = sqrtf(dst2);
                return true;
            }

            return false;
        }
    };

    


    template<typename Real>
    struct PointCloudBVH : public BVH<Real, lbvh::PointAABBGetter<Real>>
    {
        
        using ParentBVHType = BVH<Real, lbvh::PointAABBGetter<Real>>;
        using typename ParentBVHType::RealType;
        using typename ParentBVHType::AABBType;


        PointCloudBVH(bool query_host_enabled = false) :
            ParentBVHType::BVH(query_host_enabled)
        {

        }

        AabbGetterType createAabbGetter() const noexcept override
        {
			return lbvh::PointAABBGetter<Real>();
		}

        void initialize(const std::vector<Real*>& pointClouds, const std::vector<size_t>& pointCloudSizes)
        {
            assert(pointClouds.size() == pointCloudSizes.size());

            numAllPrimitives = 0;
            numGeometries = pointClouds.size();

            positionVectorsHost.resize(numGeometries);;
            primitiveVectorsHost.resize(numGeometries);;

            positionVectorsDevice.resize(numGeometries);
            primitiveVectorsDevice.resize(numGeometries);

            for (size_t iCloud = 0; iCloud < pointClouds.size(); iCloud++)
            {

                for (size_t iPt = 0; iPt < pointCloudSizes[iCloud]; iPt++) {

                    positionVectorsHost[iCloud].push_back(pointClouds[iCloud][iPt * 3]);
                    positionVectorsHost[iCloud].push_back(pointClouds[iCloud][iPt * 3 + 1]);
                    positionVectorsHost[iCloud].push_back(pointClouds[iCloud][iPt * 3 + 2]);
                    primitiveVectorsHost[iCloud].push_back(iPt);

                    primitiveInfoVectorHost.push_back(iCloud);
                    primitiveInfoVectorHost.push_back(iPt);
                }
                numAllPrimitives += pointCloudSizes[iCloud];

                // convert to device vectors
                positionVectorsDevice[iCloud] = positionVectorsHost[iCloud];
                primitiveVectorsDevice[iCloud] = primitiveVectorsHost[iCloud];

                // assemble geometry data
                geometryDataVectorHost.push_back({ 
                    (IndexType)pointCloudSizes[iCloud],
                    (IndexType)primitiveVectorsHost[iCloud].size(),
                    positionVectorsHost[iCloud].data(),
                    primitiveVectorsHost[iCloud].data()
                    });

                GeometryData<Real> geometryDataDevice = {
                    (IndexType)pointCloudSizes[iCloud],
                    (IndexType)primitiveVectorsDevice[iCloud].size(),
                    positionVectorsDevice[iCloud].data().get(),
                    primitiveVectorsDevice[iCloud].data().get()
                };

                //printf("geometryDataDevice: %d %d %p %p\n", geometryDataDevice.numPositions, geometryDataDevice.numPrimitives,
                //    geometryDataDevice.positionsBuffer, geometryDataDevice.primitivesBuffer);
                geometryDataVectorDevice.push_back(geometryDataDevice);
            }

            primitiveInfoVectorOriginalDevice = primitiveInfoVectorHost;

            geometryDataDevice = geometryDataVectorDevice.data().get();

            construct();
        };

        void setQueryPositions(const thrust::host_vector<Real>& queryPositions)
		{
            numQueryPositions = queryPositions.size() / 3;
			queryPositionsDevice = queryPositions;
            queryResultsDevice.resize(numQueryPositions);
		}

        // two points are considered to be the same if their coordinates are within a distance of 1e-6
        void queryOverlaps()
		{
            RealType* queryPositionsBuffer = queryPositionsDevice.data().get();
            GeometryData<RealType>* geometryDataBufferDeviceforLmbd = geometryDataDevice;
            IndexType* primitiveInfoBuffereOriginalDevice = primitiveInfoVectorOriginalDevice.data().get();
            const unsigned int num_internal_nodes = numAllPrimitives - 1;
            AABBType* aabbsDeviceforLmbd = aabbs_.data().get();

            thrust::for_each(thrust::device,
                thrust::make_counting_iterator<IndexType>(0),
                thrust::make_counting_iterator<IndexType>(numAllPrimitives),
                [aabbsDeviceforLmbd, 
                num_internal_nodes, 
                bvh_dev=bvh_dev,
                geometryDataBufferDeviceforLmbd]
                __device__(IndexType idx)
            {
                IndexType geomId = bvh_dev.primitiveInfoBufferDevice[idx * 2];
                IndexType primId = bvh_dev.primitiveInfoBufferDevice[idx * 2 + 1];
                Real* positionsBuffer = bvh_dev.geometryDataDevice[geomId].positionsBuffer;
                IndexType* primitivesBuffer = bvh_dev.geometryDataDevice[geomId].primitivesBuffer;
                Real* pos = positionsBuffer + 3 * *(primitivesBuffer + primId);
                printf("idx: %d, node obj_idx: %d, geomId: %d, primId: %d, nodeId: %d\npositions: %f %f %f\naabb: (%f %f %f) - (%f %f %f)\n", 
                    idx, bvh_dev.nodes[num_internal_nodes + idx].object_idx, geomId, primId, num_internal_nodes + idx,
                    pos[0], pos[1], pos[2],
                    bvh_dev.aabbs[num_internal_nodes + idx].lower.x,
                    bvh_dev.aabbs[num_internal_nodes + idx].lower.y,
                    bvh_dev.aabbs[num_internal_nodes + idx].lower.z,
                    bvh_dev.aabbs[num_internal_nodes + idx].upper.x,
                    bvh_dev.aabbs[num_internal_nodes + idx].upper.y,
                    bvh_dev.aabbs[num_internal_nodes + idx].upper.z);
            });

            thrust::for_each(thrust::device,
                thrust::make_counting_iterator<IndexType>(0),
                thrust::make_counting_iterator<IndexType>(numQueryPositions),
                [bvh_dev=bvh_dev, queryPositionsBuffer, queryResultsBuffer=queryResultsDevice.data().get()] 
                __device__(IndexType idx) {
                PointCloudClosestPtQueryResult<RealType>* queryResult = queryResultsBuffer + idx;

                queryResult->queryPosition = make_float4(queryPositionsBuffer[idx * 3],
                    queryPositionsBuffer[idx * 3 + 1], queryPositionsBuffer[idx * 3 + 2], 0.f);
                queryResult->closestPtGeomId = 0xFFFFFFFF;
                queryResult->closestPtPrimId = 0xFFFFFFFF;
                queryResult->distance = -1.f;
                const auto& self = queryResult->queryPosition;

                RealType r = 1e-6f;
                aabb<RealType> query_box;
                query_box.lower = make_float4(self.x - r, self.y - r, self.z - r, 0);
                query_box.upper = make_float4(self.x + r, self.y + r, self.z + r, 0);

                const detail::basic_device_bvh<RealType, true>& bvh = bvh_dev;

                PointOverlapCallBackFunctor<RealType> callback;

                const auto num_found = lbvh::overlap_query_device(
                    bvh, query_box, callback, queryResult);

                printf("query Id: %d, closestPtGeomId: %d, closestPtPrimId: %d, distance: %f\n", 
                    idx, queryResult->closestPtGeomId, queryResult->closestPtPrimId, queryResult->distance);

                //for (unsigned int j = 0; j < 10; ++j)
                //{
                //    const auto jdx = buffer[j];
                //    if (j >= num_found)
                //    {
                //        assert(jdx == 0xFFFFFFFF);
                //        continue;
                //    }
                //    else
                //    {
                //        assert(jdx != 0xFFFFFFFF);
                //        assert(jdx < bvh_dev.num_objects);
                //    }
                //    const auto other = bvh_dev.objects[jdx];
                //    assert(fabsf(self.x - other.x) < r); // check coordinates
                //    assert(fabsf(self.y - other.y) < r); // are in the range
                //    assert(fabsf(self.z - other.z) < r); // of query box
                //}
                //return;
            });
		}


        thrust::host_vector<thrust::host_vector<Real>> positionVectorsHost;
        thrust::host_vector<thrust::host_vector<IndexType>> primitiveVectorsHost;
        thrust::host_vector<GeometryData<Real>> geometryDataVectorHost;
        thrust::host_vector<IndexType> primitiveInfoVectorHost;

        // those should be contained in host vectors to be able to access them from the host
        thrust::host_vector<thrust::device_vector<Real>> positionVectorsDevice;
        thrust::host_vector<thrust::device_vector<IndexType>> primitiveVectorsDevice;
        thrust::device_vector<GeometryData<Real>> geometryDataVectorDevice;

        // used for closest point queries
        IndexType numQueryPositions;
        thrust::device_vector<Real> queryPositionsDevice;
        thrust::device_vector<PointCloudClosestPtQueryResult<RealType>> queryResultsDevice;

    };

}