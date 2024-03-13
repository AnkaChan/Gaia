#ifndef LBVH_BVH_CUH
#define LBVH_BVH_CUH
#include "aabb.cuh"
#include "MortonCode.h"
#include <thrust/swap.h>
#include <thrust/pair.h>
#include <thrust/tuple.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/unique.h>

namespace lbvh
{
    using IndexType = int32_t;

    template<typename Real>
    struct GeometryData {
        IndexType numPositions;
        IndexType numPrimitives;

        Real* positionsBuffer;
        IndexType* primitivesBuffer;
    };

    namespace detail
    {
        struct node
        {
            int32_t parent_idx; // parent node
            int32_t left_idx;   // index of left  child node
            int32_t right_idx;  // index of right child node
            int32_t object_idx; // == 0xFFFFFFFF if internal node.
        };

        // a set of pointers to use it on device.
        template<typename Real, bool IsConst>
        struct basic_device_bvh;
        template<typename Real>
        struct basic_device_bvh<Real, false>
        {
            using RealType = Real;
            using AABBType = aabb<RealType>;
            using NodeType = detail::node;

            unsigned int numNodes;   // (# of internal node) + (# of leaves), 2N+1
            unsigned int numObjects; // (# of leaves), the same as the number of objects

            NodeType* nodes;
            AABBType* aabbs;
            GeometryData<Real> * geometryDataDevice;
            IndexType * primitiveInfoBufferDevice; // (geometryId, primitiveId) x numPrimitives

        };
        template<typename Real>
        struct basic_device_bvh<Real, true>
        {
            using RealType = Real;
            using AABBType = aabb<RealType>;
            using NodeType = detail::node;

            unsigned int numNodes;   // (# of internal node) + (# of leaves), 2N-1
            unsigned int numObjects; // (# of leaves), the same as the number of objects

            NodeType   const* nodes;
            AABBType   const* aabbs;
            GeometryData<Real> const* geometryDataDevice;
            IndexType  const* primitiveInfoBufferDevice; // (geometryId, primitiveId) x numPrimitives
        };

        template<typename UInt>
        __device__
            inline uint2 determine_range(UInt const* node_code,
                const unsigned int num_leaves, unsigned int idx)
        {
            if (idx == 0)
            {
                return make_uint2(0, num_leaves - 1);
            }

            // determine direction of the range
            const UInt self_code = node_code[idx];
            const int L_delta = common_upper_bits(self_code, node_code[idx - 1]);
            const int R_delta = common_upper_bits(self_code, node_code[idx + 1]);
            const int d = (R_delta > L_delta) ? 1 : -1;

            // Compute upper bound for the length of the range

            const int delta_min = thrust::min(L_delta, R_delta);
            int l_max = 2;
            int delta = -1;
            int i_tmp = idx + d * l_max;
            if (0 <= i_tmp && i_tmp < num_leaves)
            {
                delta = common_upper_bits(self_code, node_code[i_tmp]);
            }
            while (delta > delta_min)
            {
                l_max <<= 1;
                i_tmp = idx + d * l_max;
                delta = -1;
                if (0 <= i_tmp && i_tmp < num_leaves)
                {
                    delta = common_upper_bits(self_code, node_code[i_tmp]);
                }
            }

            // Find the other end by binary search
            int l = 0;
            int t = l_max >> 1;
            while (t > 0)
            {
                i_tmp = idx + (l + t) * d;
                delta = -1;
                if (0 <= i_tmp && i_tmp < num_leaves)
                {
                    delta = common_upper_bits(self_code, node_code[i_tmp]);
                }
                if (delta > delta_min)
                {
                    l += t;
                }
                t >>= 1;
            }
            unsigned int jdx = idx + l * d;
            if (d < 0)
            {
                thrust::swap(idx, jdx); // make it sure that idx < jdx
            }
            return make_uint2(idx, jdx);
        }

        template<typename UInt>
        __device__
            inline unsigned int find_split(UInt const* node_code, const unsigned int num_leaves,
                const unsigned int first, const unsigned int last) noexcept
        {
            const UInt first_code = node_code[first];
            const UInt last_code = node_code[last];
            if (first_code == last_code)
            {
                return (first + last) >> 1;
            }
            const int delta_node = common_upper_bits(first_code, last_code);

            // binary search...
            int split = first;
            int stride = last - first;
            do
            {
                stride = (stride + 1) >> 1;
                const int middle = split + stride;
                if (middle < last)
                {
                    const int delta = common_upper_bits(first_code, node_code[middle]);
                    if (delta > delta_node)
                    {
                        split = middle;
                    }
                }
            } while (stride > 1);

            return split;
        }
        template<typename Real, bool IsConst, typename UInt>
        void construct_internal_nodes(const basic_device_bvh<Real, IsConst>& self,
            UInt const* node_code, const unsigned int numObjects)
        {
            thrust::for_each(thrust::device,
                thrust::make_counting_iterator<unsigned int>(0),
                thrust::make_counting_iterator<unsigned int>(numObjects - 1),
                [self, node_code, numObjects] __host__ __device__(const unsigned int idx)
            {
                self.nodes[idx].object_idx = 0xFFFFFFFF; //  internal nodes

                const uint2 ij = determine_range(node_code, numObjects, idx);
                const int gamma = find_split(node_code, numObjects, ij.x, ij.y);

                self.nodes[idx].left_idx = gamma;
                self.nodes[idx].right_idx = gamma + 1;
                if (thrust::min(ij.x, ij.y) == gamma)
                {
                    self.nodes[idx].left_idx += numObjects - 1;
                }
                if (thrust::max(ij.x, ij.y) == gamma + 1)
                {
                    self.nodes[idx].right_idx += numObjects - 1;
                }
                self.nodes[self.nodes[idx].left_idx].parent_idx = idx;
                self.nodes[self.nodes[idx].right_idx].parent_idx = idx;
                return;
            });
            return;
        }

    } // detail

    template<typename Real>
    struct default_morton_code_calculator
    {
        default_morton_code_calculator(aabb<Real> w) : whole(w) {}
        default_morton_code_calculator() = default;
        ~default_morton_code_calculator() = default;
        default_morton_code_calculator(default_morton_code_calculator const&) = default;
        default_morton_code_calculator(default_morton_code_calculator&&) = default;
        default_morton_code_calculator& operator=(default_morton_code_calculator const&) = default;
        default_morton_code_calculator& operator=(default_morton_code_calculator&&) = default;

        __device__ __host__
            inline unsigned int operator()(const aabb<Real>& box) noexcept
        {
            auto p = centroid(box);
            p.x -= whole.lower.x;
            p.y -= whole.lower.y;
            p.z -= whole.lower.z;
            p.x /= (whole.upper.x - whole.lower.x);
            p.y /= (whole.upper.y - whole.lower.y);
            p.z /= (whole.upper.z - whole.lower.z);
            return morton_code(p);
        }
        aabb<Real> whole;
    };


    template<typename Real>
    using  bvh_device = detail::basic_device_bvh<Real, false>;
    template<typename Real>
    using cbvh_device = detail::basic_device_bvh<Real, true>;
    /*
    * Real: the floating point type
    * Object
    */
    template<typename Real, typename AABBGetter,
        typename MortonCodeCalculator = default_morton_code_calculator<Real>>
        class BVH
    {
    public:
        typedef Real RealType;
        
        typedef aabb<RealType> AABBType;
        typedef detail::node NodeType;
        typedef AABBGetter AabbGetterType;
        typedef MortonCodeCalculator MortonCodeCalculatorType;

    public:

        BVH(bool query_host_enabled = false) :
            query_host_enabled_(query_host_enabled)
        {
        }


        BVH() = default;
        ~BVH() = default;
        BVH(const BVH&) = default;
        BVH(BVH&&) = default;
        BVH& operator=(const BVH&) = default;
        BVH& operator=(BVH&&) = default;

        bool  query_host_enabled() const noexcept { return query_host_enabled_; }
        bool& query_host_enabled()       noexcept { return query_host_enabled_; }

        void clear()
        {
            this->objects_h_.clear();
            this->objects_d_.clear();
            this->aabbs_h_.clear();
            this->aabbs_.clear();
            this->nodes_h_.clear();
            this->nodes_.clear();
            return;
        }

        template<typename InputIterator>
        void assign(InputIterator first, InputIterator last)
        {
            this->objects_h_.assign(first, last);
            this->objects_d_ = this->objects_h_;
            this->construct();
            return;
        }

        bvh_device<RealType> get_device_repr()       noexcept
        {
            return bvh_device<RealType>{
                static_cast<unsigned int>(nodes_.size()),
                    static_cast<unsigned int>(numAllPrimitives),
                    nodes_.data().get(), 
                    aabbs_.data().get(),
                    geometryDataDevice,
                    primitiveInfoBufferDevice

            };
        }
        cbvh_device<RealType> get_device_repr() const noexcept
        {
            return cbvh_device<RealType>{
                static_cast<unsigned int>(nodes_.size()),
                    static_cast<unsigned int>(numAllPrimitives),
                    nodes_.data().get(), 
                    aabbs_.data().get(),
                    geometryDataDevice,
                    primitiveInfoBufferDevice
            };
        }

        cbvh_device<RealType> get_device_crepr() const noexcept
        {
            return cbvh_device<RealType>{
                static_cast<unsigned int>(nodes_.size()),
                    static_cast<unsigned int>(numAllPrimitives),
                    nodes_.data().get(),
                    aabbs_.data().get(),
                    geometryDataDevice,
                    primitiveInfoBufferDevice
            };
        }

        // create your own AabbGetterType, it must be a functor, you can set the aabb construction parameters inside it
        virtual AabbGetterType createAabbGetter() const noexcept = 0;

        void construct()
        {
            if (numAllPrimitives == 0u) { return; }

            // first num_internal_nodes are none leaf nodes
            const unsigned int num_internal_nodes = numAllPrimitives - 1;
            // number of al nodes in the BVH
            const unsigned int numNodes = numAllPrimitives * 2 - 1;

            // --------------------------------------------------------------------
            // calculate morton code of each points

            const auto inf = std::numeric_limits<RealType>::infinity();
            AABBType default_aabb;
            default_aabb.upper.x = -inf; default_aabb.lower.x = inf;
            default_aabb.upper.y = -inf; default_aabb.lower.y = inf;
            default_aabb.upper.z = -inf; default_aabb.lower.z = inf;

            // resize the aabb node restorage; initialize it as default_aabb
            this->aabbs_.resize(numNodes, default_aabb);

            // this computes the aabb for each object, AabbGetterType() is a functor 
            // which takes the *InputIterator (which contains the object) and output a value
            // which is put to *OutputIterator (contains aabb)
            /*thrust::transform(this->objects_d_.begin(), this->objects_d_.end(),
                aabbs_.begin() + num_internal_nodes, AabbGetterType());*/

            AabbGetterType aabbGetter = createAabbGetter();
            GeometryData<RealType>* geometryDataBufferDeviceforLmbd = geometryDataDevice;
            IndexType* primitiveInfoBuffereOriginalDevice = primitiveInfoVectorOriginalDevice.data().get();
            AABBType* aabbsDeviceforLmbd = aabbs_.data().get();
            thrust::for_each(thrust::device,
                thrust::make_counting_iterator<IndexType>(0),
                thrust::make_counting_iterator<IndexType>(numAllPrimitives),
                [geometryDataBufferDeviceforLmbd, primitiveInfoBuffereOriginalDevice, aabbsDeviceforLmbd, num_internal_nodes, aabbGetter]
                __device__(IndexType idx)
            {
                IndexType geomId = primitiveInfoBuffereOriginalDevice[idx * 2];
                IndexType primId = primitiveInfoBuffereOriginalDevice[idx * 2 + 1];
                Real* positionsBuffer = geometryDataBufferDeviceforLmbd[geomId].positionsBuffer;
                IndexType* primitivesBuffer = geometryDataBufferDeviceforLmbd[geomId].primitivesBuffer;
                //printf("threadId: %d, geomId: %d, primId: %d, nodeId: %d\n", idx, geomId, primId, num_internal_nodes + idx);
                //printf("positionsBuffer: %p, primitivesBuffer: %p\n", positionsBuffer, primitivesBuffer);
                aabbsDeviceforLmbd[num_internal_nodes + idx] = aabbGetter(positionsBuffer, primitivesBuffer, primId);
                //aabbGetter(positionsBuffer, primitivesBuffer, idx);
            });


            // use reduced merge to get the bounding box of the whole scene, for calculating 
            // the morton codes in the next stage
            const auto aabb_whole = thrust::reduce(
                aabbs_.begin() + num_internal_nodes, aabbs_.end(), default_aabb,
                [] __host__ __device__(const AABBType & lhs, const AABBType & rhs) {
                return merge(lhs, rhs);
            });

            // caculate morton code
            thrust::device_vector<unsigned int> morton(numAllPrimitives);
            thrust::transform(aabbs_.begin() + num_internal_nodes, aabbs_.end(),
                 morton.begin(),
                MortonCodeCalculatorType(aabb_whole));

            // --------------------------------------------------------------------
            // sort object-indices by morton code

            thrust::device_vector<IndexType> indices(numAllPrimitives);
            thrust::copy(thrust::make_counting_iterator<IndexType>(0),
                thrust::make_counting_iterator<IndexType>(numAllPrimitives),
                indices.begin());
            // keep indices ascending order
            thrust::stable_sort_by_key(morton.begin(), morton.end(),
                thrust::make_zip_iterator(
                    thrust::make_tuple(aabbs_.begin() + num_internal_nodes,
                        indices.begin())));

            primitiveInfoVectorSortedDevice.resize(primitiveInfoVectorOriginalDevice.size());
            primitiveInfoBufferDevice = primitiveInfoVectorSortedDevice.data().get();
            const IndexType* indicesBufferDevice = indices.data().get();

            thrust::for_each(thrust::device,
				thrust::make_counting_iterator<IndexType>(0),
				thrust::make_counting_iterator<IndexType>(numAllPrimitives),
				[primitiveInfoBuffereOriginalDevice, indicesBufferDevice, primitiveInfoBufferDevice=primitiveInfoBufferDevice]
				__device__(IndexType idx)
			{
                primitiveInfoBufferDevice[idx * 2] = primitiveInfoBuffereOriginalDevice[indicesBufferDevice[idx] * 2];
                primitiveInfoBufferDevice[idx * 2 + 1] = primitiveInfoBuffereOriginalDevice[indicesBufferDevice[idx] * 2 + 1];
			});

            // validate the sorted primitive info
            //thrust::for_each(thrust::device,
            //    thrust::make_counting_iterator<IndexType>(0),
            //    thrust::make_counting_iterator<IndexType>(numAllPrimitives),
            //    [indicesBufferDevice, aabbsDeviceforLmbd, num_internal_nodes, geometryDataBufferDeviceforLmbd,
            //    primitiveInfoBufferDevice = primitiveInfoBufferDevice]
            //    __device__(IndexType idx)
            //{
            //    IndexType geomId = primitiveInfoBufferDevice[idx * 2];
            //    IndexType primId = primitiveInfoBufferDevice[idx * 2 + 1];
            //    Real* positionsBuffer = geometryDataBufferDeviceforLmbd[geomId].positionsBuffer;
            //    IndexType* primitivesBuffer = geometryDataBufferDeviceforLmbd[geomId].primitivesBuffer;
            //    Real* pos = positionsBuffer + 3 * *(primitivesBuffer + primId);
            //    printf("idx: %d, geomId: %d, primId: %d, nodeId: %d\npositions: %f %f %f\naabb: (%f %f %f) - (%f %f %f)\n", 
            //        idx, geomId, primId, num_internal_nodes + idx,
            //        pos[0], pos[1], pos[2],
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].lower.x,
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].lower.y,
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].lower.z,
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].upper.x,
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].upper.y,
            //        aabbsDeviceforLmbd[num_internal_nodes + idx].upper.z);
            //});

            // --------------------------------------------------------------------
            // check morton codes are unique

            thrust::device_vector<unsigned long long int> morton64(numAllPrimitives);
            const auto uniqued = thrust::unique_copy(morton.begin(), morton.end(),
                morton64.begin());

            const bool morton_code_is_unique = (morton64.end() == uniqued);
            if (!morton_code_is_unique)
            {
                thrust::transform(morton.begin(), morton.end(), indices.begin(),
                    morton64.begin(),
                    [] __host__ __device__(const unsigned int m, const unsigned int idx)
                {
                    unsigned long long int m64 = m;
                    m64 <<= 32;
                    m64 |= idx;
                    return m64;
                });
            }

            // --------------------------------------------------------------------
            // construct leaf nodes and aabbs

            NodeType default_node;
            default_node.parent_idx = 0xFFFFFFFF;
            default_node.left_idx = 0xFFFFFFFF;
            default_node.right_idx = 0xFFFFFFFF;
            default_node.object_idx = 0xFFFFFFFF;
            this->nodes_.resize(numNodes, default_node);

            //thrust::transform(indices.begin(), indices.end(),

            // primitiveInfoBufferDevice is already sorted by morton code, therefore the indices is just a sequence
            thrust::transform(thrust::make_counting_iterator<IndexType>(0),
                thrust::make_counting_iterator<IndexType>(numAllPrimitives),
                this->nodes_.begin() + num_internal_nodes,
                [] __host__ __device__(const IndexType idx)
            {
                NodeType n;
                n.parent_idx = 0xFFFFFFFF;
                n.left_idx = 0xFFFFFFFF;
                n.right_idx = 0xFFFFFFFF;
                n.object_idx = idx;
                return n;
            });

            // --------------------------------------------------------------------
            // construct internal nodes

            const auto self = this->get_device_repr();
            if (morton_code_is_unique)
            {
                const unsigned int* node_code = morton.data().get();
                detail::construct_internal_nodes(self, node_code, numAllPrimitives);
            }
            else // 64bit version
            {
                const unsigned long long int* node_code = morton64.data().get();
                detail::construct_internal_nodes(self, node_code, numAllPrimitives);
            }

            // --------------------------------------------------------------------
            // create AABB for each node by bottom-up strategy
            // this could be used for refitting the BVH

            thrust::device_vector<int> flag_container(num_internal_nodes, 0);
            const auto flags = flag_container.data().get();



            thrust::for_each(thrust::device,
                thrust::make_counting_iterator<IndexType>(num_internal_nodes),
                thrust::make_counting_iterator<IndexType>(numNodes),
                [self, flags] __device__(IndexType idx)
            {
                unsigned int parent = self.nodes[idx].parent_idx;
                while (parent != 0xFFFFFFFF) // means idx == 0
                {
                    const int old = atomicCAS(flags + parent, 0, 1);
                    if (old == 0)
                    {
                        // this is the first thread entered here.
                        // wait the other thread from the other child node.
                        return;
                    }
                    assert(old == 1);
                    // here, the flag has already been 1. it means that this
                    // thread is the 2nd thread. merge AABB of both childlen.

                    const auto lidx = self.nodes[parent].left_idx;
                    const auto ridx = self.nodes[parent].right_idx;
                    const auto lbox = self.aabbs[lidx];
                    const auto rbox = self.aabbs[ridx];
                    self.aabbs[parent] = merge(lbox, rbox);

                    // look the next parent...
                    parent = self.nodes[parent].parent_idx;
                }
                return;
            });

            bvh_dev = get_device_crepr();

            if (this->query_host_enabled_)
            {
                aabbs_h_ = aabbs_;
                nodes_h_ = nodes_;
            }


            return;
        }

        thrust::host_vector<NodeType> const& nodes_host() const noexcept { return nodes_h_; }
        thrust::host_vector<NodeType>& nodes_host()       noexcept { return nodes_h_; }
        thrust::host_vector<AABBType> const& aabbs_host() const noexcept { return aabbs_h_; }
        thrust::host_vector<AABBType>& aabbs_host()       noexcept { return aabbs_h_; }

    protected:
        thrust::host_vector  <AABBType>     aabbs_h_;
        thrust::device_vector<AABBType>     aabbs_;
        thrust::host_vector  <NodeType>     nodes_h_;
        thrust::device_vector<NodeType>     nodes_;

        // the original primitive info, which is set by the user when initializing the BVH
        thrust::device_vector<IndexType> primitiveInfoVectorOriginalDevice;


        IndexType numAllPrimitives;  // shape=(numGeometries,) number of primitives
        IndexType numGeometries;

        GeometryData<RealType>* geometryDataDevice;

        // flattened info for each primitive, each corresponds to a leaf node
        // corresponds to each leaf node
        IndexType* primitiveInfoBufferDevice; // (geometryId, primitiveId) x numPrimitives

        cbvh_device<RealType> bvh_dev;

        bool query_host_enabled_;

    private:
        // the primitive info, which is sorted by morton code, and used for BVH traversal
        // user should not access this directly
        thrust::device_vector<IndexType> primitiveInfoVectorSortedDevice;
    };

} // lbvh
#endif// LBVH_BVH_CUH
