#pragma once

#include "Eigen/core"
#include "Eigen/StdVector"
#include <array>
#include <vector>
#include <memory>
#include <mutex>

#include <MeshFrame/TetMesh/TMeshStaticLibHeaders.h>
#include <MeshFrame/TetMesh/SurfaceMesh/SurfaceMeshHeaders.h>

#include "../Types/Types.h"
#include "Materials/Materials.h"

#include "CuMatrix/Geometry/Geometry.h"
#include "CuMatrix/MatrixOps/CuMatrix.h"

// #define OUTPUT_TRAVERSED_TETS
// #define TET_TET_ADJACENT_LIST
// #define KEEP_MESHFRAME_MESHES
#define ENABLE_REST_POSE_CLOSEST_POINT

namespace GAIA {
	//using MF::TVec3Block;
	//using MF::TetMesh::Vec4BlockI;
	//using MF::IdType;
	//using MF::TetMesh::TTetIdsMat;

	typedef MF::TetMesh::CTMeshStaticDType<FloatingType> TetMeshMF;
	typedef MF::TetMesh::CSurfaceMeshStaticDType<FloatingType> TetSurfaceMeshMF;
	typedef MF::TetMesh::TIterators<TetMeshMF> TIt;

	//typedef MF::TVerticesMat<FloatingType> VerticesMat;

	enum class  TraverseStopReason
	{
		querySuccess = 0,
		overflow = 1,
		passedMaximumDis = 2,
		reachedBoundary = 3,
		emptyStack = 4

	};

	struct TraverseStatistics
	{
		int numTetsTraversed = 0;
		TraverseStopReason stopReason;
#ifdef OUTPUT_TRAVERSED_TETS 
		std::vector<int32_t>& traversedTetsOutput
#endif
	};


	inline std::string traverseStopReasonToText(TraverseStopReason& reason) {
		std::string reasonStr;
		switch (reason)
		{
		case TraverseStopReason::emptyStack:
			reasonStr = "emptyStack";
			break;
		case TraverseStopReason::overflow:
			reasonStr = "overflow";
			break;
		case TraverseStopReason::reachedBoundary:
			reasonStr = "reachedBoundary";
			break;
		case TraverseStopReason::querySuccess:
			reasonStr = "querySuccess";
			break;
		case TraverseStopReason::passedMaximumDis:
			reasonStr = "passedMaximumDis";
			break;
		default:
			break;
		}

		return reasonStr;
	}

	struct TetMeshTopology
	{

		void initialize(TetMeshMF* pTM_MF, ObjectParams::SharedPtr pObjectParams);

		typedef std::shared_ptr<TetMeshTopology> SharedPtr;
		typedef TetMeshTopology* Ptr;

		TTetIdsMat tetVIds;
		// vertex ids of the surface points
		VecDynamicI surfaceVIds;
		// -1 if not surface vertex
		VecDynamicI tetVertIndicesToSurfaceVertIndices;
		// surface faces, its vertex indices correspond to the ordering of mVerPos (not surfaceVIds)
		FaceVIdsMat surfaceFacesTetMeshVIds;
		// surface faces, its vertex indices correspond to the ordering of surfaceVIds
		FaceVIdsMat surfaceFacesSurfaceMeshVIds;
		VecDynamicI surfaceFacesBelongingTets;
		VecDynamicI surfaceFacesIdAtBelongingTets;

		// say 3 vertices in surfaceFacesTetMeshVIds are A, B and C,
		// the 3 edges will be AB, BC, CD
		// and 3 neighbor faces will be on three face on the other side of AB, BC, CA correspondingly
		FaceVIdsMat surfaceFaces3NeighborFaces;

		// - numSurfacesVerts x numFacesEachSurfaceVert
		// - - size: numSurfaceVertices; it stores surface mesh face Id
		std::vector<std::vector<IdType>> surfaceVertexNeighborSurfaceFaces;
		// - - size: numSurfaceVertices; it stores tetmesh vertex ids
		std::vector<std::vector<IdType>> surfaceVertexNeighborSurfaceVertices;

		// all vertice's neighbor tets are stacked together
		VecDynamicI vertexNeighborTets;
		// size is vertexNeighborTets.size(), stores a number ranged [0,3] for each neighbor tet, indicating which vertex of the tet is the current vertex
		VecDynamicI vertexNeighborTets_vertexOrder;
		// size is 2 * numVertices(), (2xiV)-th element stores where the neighbor tets of the info starts at
		// and (2xiV+1)-th element stores how many neighbor tets the iV-th vertex has
		VecDynamicI vertexNeighborTets_infos;

		// all edge's neighbor tets are stacked together
		VecDynamicI edgeNeighborTets;
		// size is  edgeNeighborTets.size(), stores a number ranged [0,3] for each neighbor tet, indicating the edge's 2 vertices order in corresponding 
		Mat2xI edgeNeighborTets_edge2VerticesOrderInTet;
		// size is 2 * numEdge(), (2xiV)-th element stores where the neighbor tets of the info starts at  
		// and (2xi+1)-th element stores how many neighbor tets the i-th edge has
		VecDynamicI edgeNeighborTets_infos;


#ifdef TET_TET_ADJACENT_LIST
		// for inversion solve, not just adjacient by face, but also those adjacent by faces
		std::vector<std::vector<IdType>> tetAllNeighborTets;
#endif // TET_TET_ADJACENT_LIST

		// data for collision detection & handling

		// - each tetrahedron's four vertex indices' XOR sum;
		VecDynamicI tetsXorSums;
		// - each tetrahedron's four neighbor tets, ordered by the tetrahedron cross the corresponding vertex in tetVIds;
		TTetIdsMat tetsNeighborTets;
		VecDynamicI surfaceEdges;
		VecDynamicBool tetsIsSurfaceTet;

		Mat2xI edges;
		size_t nEdges;
		size_t nVerts;

		std::vector<std::vector<int32_t>> tetsColoringCategories;
		std::vector<std::vector<int32_t>> edgesColoringCategories;
		std::vector<std::vector<int32_t>> verticesColoringCategories;

		size_t numTets() { return tetVIds.cols(); }
		size_t numVertices() { return nVerts; }
		int& getTVId(int tId, int vId) { return tetVIds(vId, tId); }

	};

	struct TetMeshFEM
	{
		typedef std::shared_ptr<TetMeshFEM> SharedPtr;
		typedef TetMeshFEM* Ptr;

		int& getTVId(int tId, int vId);

		template<typename Derived>
		void computeDs(Eigen::DenseBase<Derived>& Ds, int tId);
		Eigen::Map<Mat3> getDmInv(int tId);
		virtual void initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF);

		TVerticesMat mVertPos;
		IdType m_nVertices;
		IdType m_nTets;

		TVerticesMat mVertPrevPos;
		TVerticesMat mVelocity;

#ifdef KEEP_MESHFRAME_MESHES
		std::shared_ptr<TetMeshMF> m_pTM_MF;
#endif // KEEP_MESHFRAME_MESHES

		Eigen::Block<TVerticesMat, 3, -1> positions();
		Eigen::Block<TVerticesMat, 3, -1> positionsPrev();
		Eigen::Block<TVerticesMat, 3, -1> vertices();
		Eigen::Block<TVerticesMat, 3, -1> velocities();
		Vec3Block vertex(size_t i);
		ConstVec3Block vertex(size_t i) const;
		Vec3Block velocity(size_t i);
		Vec3Block vertexPrevPos(size_t i);
		Vec3Block surfaceVertex(size_t i);
		size_t numSurfaceVerts();
		size_t numSurfaceFaces();
		size_t numVertices();
		size_t numTets();

		void save(const std::string outFile);
		void saveAsPLY(const std::string outFile);

		IdType* getSurfaceFVIdsInTetMeshVIds(IdType surfaceFId);

		Vec4BlockI GAIA::TetMeshFEM::tet(size_t i);

		bool tetrahedralTraverseTo(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
			int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics);

		bool tetrahedralTraverseToLoopLess(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
			int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics);

		bool tetrahedralTraverseToDynamic(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
			int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics);

		int32_t getNextTet(int32_t tetId, int32_t exitFaceId);
		int exitFaceSelection(Eigen::Matrix<FloatingType, 2, 4>& ptsProj2D, Eigen::Vector4i& possibleExitFace,
			FloatingType rayTriIntersectionEpsilon);

		int TetMeshFEM::checkExitFaceForward(const Vec3& rayDir, int32_t tetId, int32_t incomingFaceIdCurTet,
			Eigen::Vector4i& possibleExitFace);

		void projectTo2DCoordinates(Eigen::Matrix<FloatingType, 2, 4>& ptsProj2D, int32_t tetId, int32_t incomingFaceId,
			const Eigen::Matrix<FloatingType, 2, 3>& axesT, const Vec3& origin);

		inline void TetMeshFEM::copyRepermutedVerts(Eigen::Matrix<FloatingType, 3, 4>& ptsPermuted3D,
			int32_t tetId, int32_t incomingFaceId);

		template <typename Derived2x1>
		inline void barycentric2DTriangle(const Eigen::MatrixBase<Derived2x1>& p, const  Eigen::MatrixBase<Derived2x1>& a,
			const  Eigen::MatrixBase<Derived2x1>& b, const  Eigen::MatrixBase<Derived2x1>& c, Vec3& barys);

		template <typename Derived2x1>
		inline void originBarycentric2DTriangle(const  Eigen::MatrixBase<Derived2x1>& a,
			const  Eigen::MatrixBase<Derived2x1>& b, const  Eigen::MatrixBase<Derived2x1>& c, Vec3& barys);

		// vId are surface vertex vId
		void computeVertexNormal(int32_t surfaceVId, Vec3& normal);
		void computeFaceOrientedVolume(int32_t surfaceFaceId, Vec3& orientedVolume);
		void computeFaceNormal(int32_t surfaceFaceId, Vec3& normal);
		// say 3 vertices in surfaceFacesTetMeshVIds[faceId] are A, B and C,
		// the 3 edges will be AB, BC, CD
		void computeEdgeNormal(int32_t faceId, int32_t edgeId, Vec3& normal);

		void computeTopology(TetMeshMF* pTM_MF);
		const VecDynamicI& surfaceVIds() const { return pTopology->surfaceVIds; }
		const FaceVIdsMat& surfaceFacesSurfaceMeshVIds() const { return pTopology->surfaceFacesSurfaceMeshVIds; }
		TTetIdsMat& tetVIds() { return pTopology->tetVIds; }
		const TTetIdsMat& tetVIds() const { return pTopology->tetVIds; }
		const VecDynamicI& tetVertIndicesToSurfaceVertIndices() { return pTopology->tetVertIndicesToSurfaceVertIndices; }
		const FaceVIdsMat& surfaceFacesTetMeshVIds() const { return pTopology->surfaceFacesTetMeshVIds; }
		const VecDynamicI& surfaceFacesBelongingTets()const { return pTopology->surfaceFacesBelongingTets; }
		const VecDynamicI& surfaceFacesIdAtBelongingTets() const { return pTopology->surfaceFacesIdAtBelongingTets; }
		const VecDynamicI& tetsXorSums() const { return pTopology->tetsXorSums; }
		const TTetIdsMat& tetsNeighborTets() const { return pTopology->tetsNeighborTets; }
		const VecDynamicI& surfaceEdges() const { return pTopology->surfaceEdges; }
		const VecDynamicBool& tetsIsSurfaceTet() const { return pTopology->tetsIsSurfaceTet; }
		const FaceVIdsMat& surfaceFaces3NeighborFaces() const { return pTopology->surfaceFaces3NeighborFaces; }
		const std::vector<std::vector<IdType>>& surfaceVertexNeighborSurfaceVertices() const { return pTopology->surfaceVertexNeighborSurfaceVertices; }
		int getNumVertexNeighborTets(int iV) const { return pTopology->vertexNeighborTets_infos(2 * iV + 1); }
		int getVertexNeighborTet(int iV, int iNeighborTet) const { return pTopology->vertexNeighborTets(pTopology->vertexNeighborTets_infos(2 * iV) + iNeighborTet); }
		int getVertexNeighborTetVertexOrder(int iV, int iNeighborTet) const { return pTopology->vertexNeighborTets_vertexOrder(pTopology->vertexNeighborTets_infos(2 * iV) + iNeighborTet); }

		size_t numEdges() { return pTopology->nEdges; }
		Mat2xI& edges() { return pTopology->edges; }
		int getNumEdgeNeighborTets(int iE) const { return pTopology->edgeNeighborTets_infos(2 * iE + 1); }
		int getEdgeNeighborTet(int iE, int iNeighborTet) const { return pTopology->edgeNeighborTets(pTopology->edgeNeighborTets_infos(2 * iE) + iNeighborTet); }
		void getEdgeNeighborTetVertexOrder(int iE, int iNeighborTet, int& corner0, int& corner1) const {
			int index = pTopology->edgeNeighborTets_infos(2 * iE) + iNeighborTet;
			corner0 = pTopology->edgeNeighborTets_edge2VerticesOrderInTet(0, index);
			corner1 = pTopology->edgeNeighborTets_edge2VerticesOrderInTet(1, index);
		}

		std::vector<std::vector<int32_t>>& tetsColoringCategories() { return pTopology->tetsColoringCategories; }
		std::vector<std::vector<int32_t>>& edgesColoringCategories() { return pTopology->edgesColoringCategories; }
		std::vector<std::vector<int32_t>>& verticesColoringCategories() { return pTopology->verticesColoringCategories; }

#ifdef TET_TET_ADJACENT_LIST
		// for inversion solve, not just adjacient by face, but also those adjacent by faces
		std::vector<std::vector<IdType>> tetAllNeighborTets;
#endif // TET_TET_ADJACENT_LIST

		// nTets x 9: restpose DsInvs flattened in col major and concatenated
		VecDynamic DmInvs;

		// per-object parameters
		ObjectParams::SharedPtr pObjectParams;

		bool DCDEnabled(int32_t iV);
		bool CCDEnabled(int32_t iV);

		//void computeIntersectionPoint(const Vec3& barys, int32_t currentTetId, int32_t incomingFaceId,
		//	int32_t exitFaceId, Vec3& intersectionPt);

		VecDynamic tetRestVolume;
		VecDynamic tetInvRestVolume;
		VecDynamic vertexMass;
		VecDynamic vertexInvMass;

		// vertex ids are tetmesh ids
		VecDynamicBool tetsInvertedSign;
		VecDynamicBool verticesInvertedSign;

		VecDynamicBool verticesInvertedSignPrevPos;
		VecDynamicBool tetsInvertedSignPrevPos;

		// vertex ids are tetmesh ids
		VecDynamicBool verticesCollisionDetectionEnabled;

		static const int32_t passThroughCheckSteps;
		static const int32_t tet4Faces[4][3];
		// how the 3 possible exit face is made of; note that the indices corresponds to ptsProj2D 
		// (which have been re-permuted such that the incoming face will be f3
		static const int32_t posibleExit3Faces[3][3];

#ifdef ENABLE_REST_POSE_CLOSEST_POINT
		// restpose infos
		TVerticesMat restposeVerts;
#endif // DEBUG

		VecDynamicBool fixedMask;

		// those are for CCD, to apply CDD the vertex must be inversion free in both prev position and current position
		bool activeForCollision = true;
		bool activeForMaterialSolve = true;

		// its id in the corresponding physic framework
		IdType meshId = -1;

		TetMeshTopology::SharedPtr pTopology;
		// records all the topology of the previously loaded mesh
		// we can grab reused the one that's already been computed since they don't change over time
		static std::map<std::string, TetMeshTopology::SharedPtr> topologies;
		static std::mutex topologies_lock;
	};


	template<typename Derived>
	inline void TetMeshFEM::computeDs(Eigen::DenseBase<Derived>& Ds, int tId)
	{
		auto v1 = mVertPos.block<3, 1>(0, tetVIds()(0, tId));
		for (size_t iCol = 0; iCol < 3; iCol++)
		{
			Ds.block<3, 1>(0, iCol) = mVertPos.block<3, 1>(0, tetVIds()(iCol + 1, tId)) - v1;
		}
	}

	inline Eigen::Block<TVerticesMat, 3, -1> TetMeshFEM::vertices()
	{
		return mVertPos.block<3, -1>(0, 0, 3, numVertices());;
	}

	inline Eigen::Block<TVerticesMat, 3, -1> TetMeshFEM::positions()
	{
		return mVertPos.block<3, -1>(0, 0, 3, numVertices());;
	}

	inline Eigen::Block<TVerticesMat, 3, -1> TetMeshFEM::positionsPrev()
	{
		return mVertPrevPos.block<3, -1>(0, 0, 3, numVertices());;
	}

	inline Eigen::Block<TVerticesMat, 3, -1> TetMeshFEM::velocities()
	{
		return mVelocity.block<3, -1>(0, 0, 3, numVertices());;
	}

	inline Vec3Block GAIA::TetMeshFEM::vertex(size_t i)
	{
		return mVertPos.block<3, 1>(0, i);
	}

	inline ConstVec3Block GAIA::TetMeshFEM::vertex(size_t i) const
	{
		return mVertPos.block<3, 1>(0, i);
	}


	inline Vec3Block GAIA::TetMeshFEM::velocity(size_t i)
	{
		return mVelocity.block<3, 1>(0, i);
	}

	inline Vec3Block GAIA::TetMeshFEM::vertexPrevPos(size_t i) {
		return mVertPrevPos.block<3, 1>(0, i);

	};


	inline Vec3Block GAIA::TetMeshFEM::surfaceVertex(size_t i)
	{
		return vertex(pTopology->surfaceVIds(i));
	}


	inline Vec4BlockI GAIA::TetMeshFEM::tet(size_t i)
	{
		return pTopology->tetVIds.block<4, 1>(0, i);
	}

	inline size_t GAIA::TetMeshFEM::numSurfaceVerts() {
		return pTopology->surfaceVIds.size();
	};

	inline size_t GAIA::TetMeshFEM::numSurfaceFaces() {
		return pTopology->surfaceFacesSurfaceMeshVIds.cols();
	}

	inline int& GAIA::TetMeshFEM::getTVId(int tId, int vId)
	{
		return pTopology->tetVIds(vId, tId);
	}


	inline Eigen::Map<Mat3> GAIA::TetMeshFEM::getDmInv(int tId)
	{
		return Eigen::Map<Mat3>(DmInvs.data() + tId * 9);
	}

	inline int32_t GAIA::TetMeshFEM::getNextTet(int32_t tetId, int32_t exitFaceId)
	{
		return int32_t();
	}

	inline FloatingType signedSquare(FloatingType v) { return copysignf(v * v, v); }


	inline int TetMeshFEM::exitFaceSelection(Eigen::Matrix<FloatingType, 2, 4>& ptsProj2D,
		Eigen::Vector4i& possibleExitFace, FloatingType rayTriIntersectionEpsilon)
	{
		int numExitFaces = 0;
		//rayTriIntersectionEpsilon = -rayTriIntersectionEpsilon;
		possibleExitFace = Eigen::Vector4i::Zero();

		Eigen::Matrix<FloatingType, 2, 1> v1 = ptsProj2D.col(1) - ptsProj2D.col(0);
		Eigen::Matrix<FloatingType, 2, 1> v2 = ptsProj2D.col(2) - ptsProj2D.col(0);

		FloatingType inComingTriangleArea = CuMatrix::vec2CrossProduct(v1.data(), v2.data());
		FloatingType inComingTriangleAreaSign = copysignf(1.0f, inComingTriangleArea);

		//FloatingType detP3P0Square = signedSquare(CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(0).data()));
		//FloatingType detP3P1Square = signedSquare(CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(1).data()));
		//FloatingType detP3P2Square = signedSquare(CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(2).data()));

		FloatingType detP3P0 = CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(0).data());
		FloatingType detP3P1 = CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(1).data());
		FloatingType detP3P2 = CuMatrix::vec2CrossProduct(ptsProj2D.col(3).data(), ptsProj2D.col(2).data());


		if (detP3P1 * inComingTriangleAreaSign >= -rayTriIntersectionEpsilon
			&& detP3P2 * inComingTriangleAreaSign <= rayTriIntersectionEpsilon)
		{
			possibleExitFace(0) = 1;
			++numExitFaces;
		}

		if (detP3P2 * inComingTriangleAreaSign >= -rayTriIntersectionEpsilon
			&& detP3P0 * inComingTriangleAreaSign <= rayTriIntersectionEpsilon)
		{
			possibleExitFace(1) = 1;
			++numExitFaces;
		}

		if (detP3P0 * inComingTriangleAreaSign >= -rayTriIntersectionEpsilon
			&& detP3P1 * inComingTriangleAreaSign <= rayTriIntersectionEpsilon)
		{
			possibleExitFace(2) = 1;
			++numExitFaces;
		}

		//FloatingType p0Norm2 = ptsProj2D.col(0).squaredNorm();
		//FloatingType p1Norm2 = ptsProj2D.col(1).squaredNorm();
		//FloatingType p2Norm2 = ptsProj2D.col(2).squaredNorm();
		//FloatingType p3Norm2 = ptsProj2D.col(3).squaredNorm();
		//if (detP3P1 * inComingTriangleAreaSign >= -(p3Norm2 + p1Norm2) * rayTriIntersectionEpsilon
		//	&& detP3P2 * inComingTriangleAreaSign <= (p3Norm2 + p2Norm2) * rayTriIntersectionEpsilon)
		//{
		//	possibleExitFace(0) = 1;
		//}

		//if (detP3P2 * inComingTriangleAreaSign >= -(p3Norm2 + p2Norm2) * rayTriIntersectionEpsilon
		//	&& detP3P0 * inComingTriangleAreaSign <= (p3Norm2 + p0Norm2) * rayTriIntersectionEpsilon)
		//{
		//	possibleExitFace(1) = 1;
		//}

		//if (detP3P0 * inComingTriangleAreaSign >= -(p3Norm2 + p0Norm2) * rayTriIntersectionEpsilon
		//	&& detP3P1 * inComingTriangleAreaSign <= (p3Norm2 + p1Norm2) * rayTriIntersectionEpsilon)
		//{
		//	possibleExitFace(2) = 1;
		//}

		return numExitFaces;
	}

	inline int TetMeshFEM::checkExitFaceForward(const Vec3& rayDir, int32_t currentTetId, int32_t incomingFaceIdCurTet,
		Eigen::Vector4i& possibleExitFace)
	{
		Eigen::Matrix<FloatingType, 3, 4> ptsPermuted3D;
		copyRepermutedVerts(ptsPermuted3D, currentTetId, incomingFaceIdCurTet);
		int32_t* tetVIds = pTopology->tetVIds.col(currentTetId).data();

		float tetOrientedVolume = CuMatrix::tetOrientedVolume(mVertPos.data(), tetVIds);
		FloatingType tetOrientedVolumeSign = copysignf(1.0f, tetOrientedVolume);

		int numExitFaces = 0;
		for (size_t iF = 0; iF < 3; iF++)
		{
			if (possibleExitFace(iF)) {
				Vec3 exitFaceNormal;
				int32_t exitFaceId = tet4Faces[incomingFaceIdCurTet][iF];
				CuMatrix::triangleOrientedArea(mVertPos.data(), tetVIds[tet4Faces[incomingFaceIdCurTet][0]],
					tetVIds[tet4Faces[incomingFaceIdCurTet][1]], tetVIds[tet4Faces[incomingFaceIdCurTet][2]], exitFaceNormal.data());

				if (rayDir.dot(exitFaceNormal / exitFaceNormal.norm()) * tetOrientedVolumeSign > 0)
				{
					//possibleExitFace(iF) = 0;
				}
				else
				{
					numExitFaces++;
				}
			}
		}

		return numExitFaces;
	}


	inline void TetMeshFEM::projectTo2DCoordinates(Eigen::Matrix<FloatingType, 2, 4>& ptsProj2D, int32_t tetId, int32_t incomingFaceId,
		const Eigen::Matrix<FloatingType, 2, 3>& axesT, const Vec3& origin)
	{
		int32_t* tetVIds = pTopology->tetVIds.col(tetId).data();
		int32_t repermutedTetVIds[4] = {
			tetVIds[tet4Faces[incomingFaceId][0]],
			tetVIds[tet4Faces[incomingFaceId][1]],
			tetVIds[tet4Faces[incomingFaceId][2]],
			tetVIds[incomingFaceId]
		};

		for (int iV = 0; iV < 4; iV++)
		{
			// ptsPermuted3D.col(iV) = vertex(repermutedTetVIds[iV]);
			ptsProj2D.col(iV) = axesT * (vertex(repermutedTetVIds[iV]) - origin);
		}

	}

	inline void TetMeshFEM::copyRepermutedVerts(Eigen::Matrix<FloatingType, 3, 4>& ptsPermuted3D,
		int32_t tetId, int32_t incomingFaceId)
	{
		int32_t* tetVIds = pTopology->tetVIds.col(tetId).data();
		int32_t repermutedTetVIds[4] = {
			tetVIds[tet4Faces[incomingFaceId][0]],
			tetVIds[tet4Faces[incomingFaceId][1]],
			tetVIds[tet4Faces[incomingFaceId][2]],
			tetVIds[incomingFaceId]
		};

		for (int iV = 0; iV < 4; iV++)
		{
			ptsPermuted3D.col(iV) = vertex(repermutedTetVIds[iV]);
		}

	}

	inline void TetMeshFEM::computeFaceOrientedVolume(int32_t surfaceFaceId, Vec3& orientedVolume)
	{
		Vec3 AB = mVertPos.col(pTopology->surfaceFacesTetMeshVIds(1, surfaceFaceId)) - mVertPos.col(pTopology->surfaceFacesTetMeshVIds(0, surfaceFaceId));
		Vec3 AC = mVertPos.col(pTopology->surfaceFacesTetMeshVIds(2, surfaceFaceId)) - mVertPos.col(pTopology->surfaceFacesTetMeshVIds(0, surfaceFaceId));

		orientedVolume = AB.cross(AC);
	}

	inline void TetMeshFEM::computeFaceNormal(int32_t surfaceFaceId, Vec3& normal)
	{
		computeFaceOrientedVolume(surfaceFaceId, normal);
		normal /= normal.norm();

	}

	inline void TetMeshFEM::computeEdgeNormal(int32_t faceId, int32_t edgeId, Vec3& normal)
	{
		int32_t neiFaceId = pTopology->surfaceFaces3NeighborFaces(edgeId, faceId);

		Vec3 n1, n2;

		computeFaceOrientedVolume(faceId, n1);
		computeFaceOrientedVolume(neiFaceId, n2);
		normal = n1 + n2;
		normal /= normal.norm();
	}

	inline bool TetMeshFEM::DCDEnabled(int32_t iV)
	{
		return !verticesInvertedSign[iV] && verticesCollisionDetectionEnabled[iV];
	}

	inline bool TetMeshFEM::CCDEnabled(int32_t iV)
	{
		return !verticesInvertedSign[iV] && verticesCollisionDetectionEnabled[iV];
	}

	inline void TetMeshFEM::computeVertexNormal(int32_t surfaceVId, Vec3& normal)
	{
		normal = Vec3::Zero();
		for (int iF = 0; iF < pTopology->surfaceVertexNeighborSurfaceFaces[surfaceVId].size(); iF++)
		{
			// compute face area
			Vec3 faceOrientedVolume;
			computeFaceOrientedVolume(pTopology->surfaceVertexNeighborSurfaceFaces[surfaceVId][iF], faceOrientedVolume);

			normal += faceOrientedVolume;
		}

		normal /= normal.norm();
	}

	template <typename Derived2x1>
	inline void TetMeshFEM::barycentric2DTriangle(const Eigen::MatrixBase<Derived2x1>& p, const  Eigen::MatrixBase<Derived2x1>& a,
		const  Eigen::MatrixBase<Derived2x1>& b, const  Eigen::MatrixBase<Derived2x1>& c, Vec3& barys)
	{
		Vec2 v0 = b - a, v1 = c - a, v2 = p - a;
		FloatingType den = v0[0] * v1[1] - v1[0] * v0[1];
		barys[1] = (v2[0] * v1[1] - v1[0] * v2[1]) / den;
		barys[2] = (v0[0] * v2[1] - v2[0] * v0[1]) / den;
		barys[0] = 1.0f - barys[1] - barys[2];
	}

	template<typename Derived2x1>
	inline void TetMeshFEM::originBarycentric2DTriangle(const Eigen::MatrixBase<Derived2x1>& a, const Eigen::MatrixBase<Derived2x1>& b, const Eigen::MatrixBase<Derived2x1>& c, Vec3& barys)
	{
		Vec2 v0 = b - a, v1 = c - a, v2 = -a;
		FloatingType den = v0[0] * v1[1] - v1[0] * v0[1];
		barys[1] = (v2[0] * v1[1] - v1[0] * v2[1]) / den;
		barys[2] = (v0[0] * v2[1] - v2[0] * v0[1]) / den;
		barys[0] = 1.0f - barys[1] - barys[2];
	}

	inline IdType* GAIA::TetMeshFEM::getSurfaceFVIdsInTetMeshVIds(IdType surfaceFId)
	{
		return pTopology->surfaceFacesTetMeshVIds.col(surfaceFId).data();
	}
}