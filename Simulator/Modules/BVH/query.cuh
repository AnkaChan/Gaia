#ifndef LBVH_QUERY_CUH
#define LBVH_QUERY_CUH
#include "predicator.cuh"
#include <CuMatrix/MatrixOps/CuMatrix.h>

namespace lbvh
{
    template<typename Real>
    struct QueryDataBase
    {
        typename vector_of_t<Real> queryPosition;

    };


// query object indices that potentially overlaps with query aabb.
//
// requirements:
// - Callback a callable object that takes (obj_idx, geomId, primId, bvh, userData) as arguments, 
//   and return a boolean value that indicates whether the intersection exists or not.
// - Input is a bounding box, whose size cannot change.
//
template<typename Real, typename Callback>
__device__
unsigned int overlap_query_device(
    const detail::basic_device_bvh<Real, true>& bvh,
    const aabb<Real>& q, const Callback & overlapCallBack, QueryDataBase<Real>* queryData) noexcept
{
    using bvh_type   = detail::basic_device_bvh<Real, true>;
    using AABBType  = typename bvh_type::AABBType;
    using NodeType  = typename bvh_type::NodeType;

    IndexType  stack[64]; // is it okay?
    IndexType* stack_ptr = stack;
    *stack_ptr++ = 0; // root node is always 0

    unsigned int num_found = 0;
    do
    {
        const IndexType node  = *--stack_ptr;
        const IndexType L_idx = bvh.nodes[node].left_idx;
        const IndexType R_idx = bvh.nodes[node].right_idx;

        if(intersects(q, bvh.aabbs[L_idx]))
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                //*outiter++ = obj_idx;
                IndexType geomId = bvh.primitiveInfoBufferDevice[obj_idx * 2];
                IndexType primId = bvh.primitiveInfoBufferDevice[obj_idx * 2 + 1];
                if (overlapCallBack(obj_idx, geomId, primId, bvh, queryData))
				{
					++num_found;
				}
            }
            else if (stack_ptr < stack + 63)
            // the node is not a leaf.
            {
                *stack_ptr++ = L_idx;
            }
        }
        if(intersects(q, bvh.aabbs[R_idx]))
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                //*outiter++ = obj_idx;
                IndexType geomId = bvh.primitiveInfoBufferDevice[obj_idx * 2];
                IndexType primId = bvh.primitiveInfoBufferDevice[obj_idx * 2 + 1];
                if (overlapCallBack(obj_idx, geomId, primId, bvh, queryData))
                {
                    ++num_found;
                }
            }
            else if (stack_ptr < stack + 63)
            // the node is not a leaf.
            {
                *stack_ptr++ = R_idx;
            }
        }
    }
    while (stack < stack_ptr);
    return num_found;
}

// query object index that is the nearst to the query point.
//
// requirements:
// - DistanceCalculator must be able to calc distance between a point to an object.
//
template<typename Real, bool IsConst,
         typename DistanceCalculator>
__device__
thrust::pair<unsigned int, Real> query_device(
        const detail::basic_device_bvh<Real, IsConst>& bvh,
        const query_nearest<Real>& q, DistanceCalculator calc_dist) noexcept
{
    using bvh_type   = detail::basic_device_bvh<Real, IsConst>;
    using RealType  = typename bvh_type::RealType;
    using index_type = typename bvh_type::index_type;
    using AABBType  = typename bvh_type::AABBType;
    using NodeType  = typename bvh_type::NodeType;

    // pair of {node_idx, mindist}
    thrust::pair<index_type, RealType>  stack[64];
    thrust::pair<index_type, RealType>* stack_ptr = stack;
    *stack_ptr++ = thrust::make_pair(0, mindist(bvh.aabbs[0], q.target));

    unsigned int nearest = 0xFFFFFFFF;
    RealType dist_to_nearest_object = infinity<RealType>();
    do
    {
        const auto node  = *--stack_ptr;
        if(node.second > dist_to_nearest_object)
        {
            // if aabb mindist > already_found_mindist, it cannot have a nearest
            continue;
        }

        const index_type L_idx = bvh.nodes[node.first].left_idx;
        const index_type R_idx = bvh.nodes[node.first].right_idx;

        const AABBType& L_box = bvh.aabbs[L_idx];
        const AABBType& R_box = bvh.aabbs[R_idx];

        const RealType L_mindist = mindist(L_box, q.target);
        const RealType R_mindist = mindist(R_box, q.target);

        const RealType L_minmaxdist = minmaxdist(L_box, q.target);
        const RealType R_minmaxdist = minmaxdist(R_box, q.target);

        // there should be an object that locates within minmaxdist.

        if(L_mindist <= R_minmaxdist) // L is worth considering
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const RealType dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if(dist <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                *stack_ptr++ = thrust::make_pair(L_idx, L_mindist);
            }
        }
        if(R_mindist <= L_minmaxdist) // R is worth considering
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const RealType dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if(dist <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                *stack_ptr++ = thrust::make_pair(R_idx, R_mindist);
            }
        }
        assert(stack_ptr < stack + 64);
    }
    while (stack < stack_ptr);
    return thrust::make_pair(nearest, dist_to_nearest_object);
}
// query object index that is within a certain radius.
//
// requirements:
// - DistanceCalculator must be able to calc distance between a point to an object.
//
template<typename Real, bool IsConst,
    typename DistanceCalculator, typename OutputIterator>
__device__
unsigned int radius_query_device(
    const detail::basic_device_bvh<Real, IsConst>& bvh,
    const query_nearest<Real>& q, int32_t selfId, DistanceCalculator calc_dist,
    Real r, OutputIterator outiter,
    unsigned int max_buffer_size) noexcept
{
    using bvh_type = detail::basic_device_bvh<Real, IsConst>;
    using RealType = typename bvh_type::RealType;
    using index_type = typename bvh_type::index_type;
    using AABBType = typename bvh_type::AABBType;
    using NodeType = typename bvh_type::NodeType;

    // pair of {node_idx, mindist}
    index_type  stack[64]; // is it okay?
    index_type* stack_ptr = stack;
    *stack_ptr++ = 0;

    unsigned int num_found = 0;

    do
    {
        const auto node = *--stack_ptr;
        const index_type L_idx = bvh.nodes[node].left_idx;
        const index_type R_idx = bvh.nodes[node].right_idx;

        const AABBType& L_box = bvh.aabbs[L_idx];
        const AABBType& R_box = bvh.aabbs[R_idx];

        const RealType L_mindist = mindist(L_box, q.target);
        const RealType R_mindist = mindist(R_box, q.target);

        // there should be an object that locates within minmaxdist.

        if (L_mindist <= r) // L is worth considering
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            
            if (obj_idx != 0xFFFFFFFF) // leaf node
            {
                if (obj_idx == selfId)
                {
                    continue;
                }
                const RealType dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if (dist <= r)
                {
                    *outiter++ = obj_idx;
                    ++num_found;
                }

                if (num_found >= max_buffer_size)
                {
                    break;
                }
            }
            else
            {
                *stack_ptr++ = L_idx;
            }
        }
        if (R_mindist <= r) // R is worth considering
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if (obj_idx != 0xFFFFFFFF) // leaf node
            {
                if (obj_idx == selfId)
                {
                    continue;
                }
                const RealType dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if (dist <= r)
                {
                    *outiter++ = obj_idx;
                    ++num_found;
                }

                if (num_found >= max_buffer_size)
                {
                    break;
                }
            }
            else
            {
                *stack_ptr++ = R_idx;
            }
        }
        assert(stack_ptr < stack + 64);
    } while (stack < stack_ptr);
    return num_found;
}
//
//template<typename Real, typename Objects, typename AABBGetter,
//         typename MortonCodeCalculator, typename OutputIterator>
//unsigned int query_host(
//    const BVH<Real, Objects, AABBGetter, MortonCodeCalculator>& tree,
//    const query_overlap<Real> q, OutputIterator outiter,
//    const unsigned int max_buffer_size = 0xFFFFFFFF)
//{
//    using bvh_type   = ::lbvh::bvh<Real, Objects, AABBGetter, MortonCodeCalculator>;
//    using index_type = typename bvh_type::index_type;
//    using AABBType  = typename bvh_type::AABBType;
//    using NodeType  = typename bvh_type::NodeType;
//
//    if(!tree.query_host_enabled())
//    {
//        throw std::runtime_error("lbvh::bvh query_host is not enabled");
//    }
//
//    std::vector<std::size_t> stack;
//    stack.reserve(64);
//    stack.push_back(0);
//
//    unsigned int num_found = 0;
//    do
//    {
//        const index_type node  = stack.back(); stack.pop_back();
//        const index_type L_idx = tree.nodes_host()[node].left_idx;
//        const index_type R_idx = tree.nodes_host()[node].right_idx;
//
//        if(intersects(q.target, tree.aabbs_host()[L_idx]))
//        {
//            const auto obj_idx = tree.nodes_host()[L_idx].object_idx;
//            if(obj_idx != 0xFFFFFFFF)
//            {
//                if(num_found < max_buffer_size)
//                {
//                    *outiter++ = obj_idx;
//                }
//                ++num_found;
//            }
//            else // the node is not a leaf.
//            {
//                stack.push_back(L_idx);
//            }
//        }
//        if(intersects(q.target, tree.aabbs_host()[R_idx]))
//        {
//            const auto obj_idx = tree.nodes_host()[R_idx].object_idx;
//            if(obj_idx != 0xFFFFFFFF)
//            {
//                if(num_found < max_buffer_size)
//                {
//                    *outiter++ = obj_idx;
//                }
//                ++num_found;
//            }
//            else // the node is not a leaf.
//            {
//                stack.push_back(R_idx);
//            }
//        }
//    }
//    while (!stack.empty());
//    return num_found;
//}
//
//template<typename Real, typename AABBGetter,
//         typename MortonCodeCalculator, typename DistanceCalculator>
//std::pair<unsigned int, Real> query_host(
//        const BVH<Real, AABBGetter, MortonCodeCalculator>& tree,
//        const query_nearest<Real>& q, DistanceCalculator calc_dist) noexcept
//{
//    using bvh_type   = ::lbvh::bvh<Real, Objects, AABBGetter, MortonCodeCalculator>;
//    using RealType  = typename bvh_type::RealType;
//    using index_type = typename bvh_type::index_type;
//    using AABBType  = typename bvh_type::AABBType;
//    using NodeType  = typename bvh_type::NodeType;
//
//    if(!tree.query_host_enabled())
//    {
//        throw std::runtime_error("lbvh::bvh query_host is not enabled");
//    }
//
//    // pair of {node_idx, mindist}
//    std::vector<std::pair<index_type, RealType>> stack = {
//        {0, mindist(tree.aabbs_host()[0], q.target)}
//    };
//    stack.reserve(64);
//
//    unsigned int nearest = 0xFFFFFFFF;
//    RealType current_nearest_dist = infinity<RealType>();
//    do
//    {
//        const auto node = stack.back(); stack.pop_back();
//        if(node.second > current_nearest_dist)
//        {
//            // if aabb mindist > already_found_mindist, it cannot have a nearest
//            continue;
//        }
//
//        const index_type L_idx = tree.nodes_host()[node.first].left_idx;
//        const index_type R_idx = tree.nodes_host()[node.first].right_idx;
//
//        const AABBType& L_box = tree.aabbs_host()[L_idx];
//        const AABBType& R_box = tree.aabbs_host()[R_idx];
//
//        const RealType L_mindist = mindist(L_box, q.target);
//        const RealType R_mindist = mindist(R_box, q.target);
//
//        const RealType L_minmaxdist = minmaxdist(L_box, q.target);
//        const RealType R_minmaxdist = minmaxdist(R_box, q.target);
//
//       // there should be an object that locates within minmaxdist.
//
//        if(L_mindist <= R_minmaxdist) // L is worth considering
//        {
//            const auto obj_idx = tree.nodes_host()[L_idx].object_idx;
//            if(obj_idx != 0xFFFFFFFF) // leaf node
//            {
//                const RealType dist = calc_dist(q.target, tree.objects_host()[obj_idx]);
//                if(dist <= current_nearest_dist)
//                {
//                    current_nearest_dist = dist;
//                    nearest = obj_idx;
//                }
//            }
//            else
//            {
//                stack.emplace_back(L_idx, L_mindist);
//            }
//        }
//        if(R_mindist <= L_minmaxdist) // R is worth considering
//        {
//            const auto obj_idx = tree.nodes_host()[R_idx].object_idx;
//            if(obj_idx != 0xFFFFFFFF) // leaf node
//            {
//                const RealType dist = calc_dist(q.target, tree.objects_host()[obj_idx]);
//                if(dist <= current_nearest_dist)
//                {
//                    current_nearest_dist = dist;
//                    nearest = obj_idx;
//                }
//            }
//            else
//            {
//                stack.emplace_back(R_idx, R_mindist);
//            }
//        }
//    }
//    while (!stack.empty());
//    return std::make_pair(nearest, current_nearest_dist);
//}
} // lbvh
#endif// LBVH_QUERY_CUH
