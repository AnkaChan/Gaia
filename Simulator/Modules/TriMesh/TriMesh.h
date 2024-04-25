#pragma once

#include "Eigen/core"
#include "Eigen/StdVector"
#include <array>
#include <vector>
#include <memory>
#include <mutex>
#include <map>
#include <string>

#include "../Types/Types.h"
#include "Materials/Materials.h"

#include "CuMatrix/Geometry/Geometry.h"
#include "CuMatrix/MatrixOps/CuMatrix.h"
#include <MeshFrame/Memory/Array.h>


#define TRIANGLE_INTERNAL_FORCE(mat, iV, iFace) (mat.block<3,1>(iV * 3, iFace))
#define EDGE_INTERNAL_FORCE(mat, iV, iEdge) (mat.block<3,1>(iV * 3, iEdge))

namespace GAIA {
	struct TriMeshFEM;

	struct TriMeshParams : public ObjectParams {
		typedef std::shared_ptr<TriMeshParams> SharedPtr;
		typedef TriMeshParams* Ptr;

		std::string triangleColoringCategoriesPath;
		std::string initialState;

		bool use3DRestpose = true;

		FloatingType UVScale = 150.f;

		FloatingType frictionDynamic = 0.1f;
		FloatingType frictionEpsV = 0.01f;

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};

	struct EdgeInfo {
		int eV1, eV2,
			eV12Next, // the vertex opposite to edge (v1, v2), in the nei triangle whose orientation matches that order
			eV21Next; // the vertex opposite to edge (v2, v1), in the nei triangle whose orientation matches that order

		int fId1,  //  the nei triangle whose orientation matches (v1, v2)
			fId2;  //  the nei triangle whose orientation matches (v2, v1)

		int eV1FaceOrder[2], eV2FaceOrder[2]; // the order of the edge in the face's vertex list
	};

	struct TriMeshTopology
	{

		virtual void initialize(TriMeshFEM* pTriMesh, TriMeshParams* pObjectParams);

		typedef std::shared_ptr<TriMeshTopology> SharedPtr;
		typedef TriMeshTopology* Ptr;

		size_t numEdges = 0;
		std::vector<EdgeInfo> edgeInfos;

		// say 3 vertices in face are A, B and C,
		// the 3 edges will be AB, BC, CA
		// and 3 neighbor faces will be on three face on the other side of AB, BC, CA correspondingly
		FaceVIdsMat faces3NeighborEdges;
		FaceVIdsMat faces3NeighborFaces;

		// all vertice's neighbor faces are stacked together
		VecDynamicI vertexNeighborFaces;
		// size is numVertices(), range [0,2], stores the vertex's order in the face's vertex list 
		VecDynamicI vertexNeighborFaces_vertexOrder;
		// size is 2 x numVertices(), (2xiV)-th element stores where the neighbor faces of the info starts at  faces3NeighborFaces
		// and (2xiV+1)-th element stores how many neighbor faces the iV-th vertex has
		VecDynamicI vertexNeighborFaces_infos;

		// all vertice's neighbor edges are stacked together
		VecDynamicI vertexNeighborEdges;
		// size is numVertices(), range[0, 1], stores the vertex's order in the face's vertex list
		VecDynamicI vertexNeighborEdges_vertexOrder;
		// size is 2 x numVertices(), (2xiV)-th element stores where the neighbor edges of the info starts at vertexNeighborEdges
		// and (2xiV+1)-th element stores how many neighbor faces the iV-th vertex has
		VecDynamicI vertexNeighborEdges_infos;

		// all vertice's relavant bending are stacked together
		// each bending energy relates to 4 vertices, as opposite to 2 vertices on the edge
		VecDynamicI vertexRelevantBendings;
		// size is numVertices(), range [0,3], stores the vertex's order in the bending energy, the number corresponds to
		// 0: eV1, 1: eV2, 2: eV12Next, 3: eV21Next
		VecDynamicI vertexRelevantBendings_vertexOrder;
		// size is 2 x numVertices(), (2xiV)-th element stores where the relavant bending energy of the info starts at vertexRelavantBendings
		// and (2xiV+1)-th element stores how many neighbor faces the iV-th vertex has
		VecDynamicI vertexRelevantBendings_infos;

		// all vertice's neighbor surface vertices are stacked together
		VecDynamicI vertexNeighborVertices;
		// similiar to vertexNeighborFaces_infos
		VecDynamicI vertexNeighborVertices_infos;
		std::vector<std::vector<int32_t>> verticesColoringCategories;
	};

	struct TriMeshFEM
	{
		typedef std::shared_ptr<TriMeshFEM> SharedPtr;
		typedef TriMeshFEM* Ptr;

		void initialize(TriMeshParams::SharedPtr inObjectParams, bool precomputeToplogy);
		void applyRotationScalingTranslation(); // in thus an order
		void computeTopology();
		// if you need to use customized topology information, inherit TriMeshTopology create your own topology class
		// then override this function to return the pointer to your own topology class
		virtual TriMeshTopology::SharedPtr createTopology();

		void loadObj(std::string objFile);
		void saveAsPLY(std::string objFile);

		Eigen::Block<TVerticesMat, 3, -1> positions();

		int numFaces();
		int numEdges();
		int numVertices();
		int facePosVId(int tId, int vId) const;
		void computeDsFromPosition(int fId, Mat3x2& Ds);
		void calculateDeformationGradient(int iF, Mat3x2& F);
		void computeRestposeTrianglesFrom3D();
		void computeRestposeTrianglesFrom2D();
		Eigen::Map<Mat2> getDmInv(int fId);

		Vec3 computeNormal(int fId, bool normalize = true) const;

		void evaluateFaceNormals(bool normalize = true);

		void get2DProjectionAxes(Vec3& normal, Eigen::Matrix<FloatingType, 3, 2>& axesT);
		void projectTriangleTo2D(int iF, const Eigen::Matrix<FloatingType, 3, 2>& axesT, Eigen::Matrix<FloatingType, 2, 3>& tri2D);
		Vec3Block vertex(size_t i);
		ConstVec3Block vertex(size_t i) const;
		ConstVec3Block vertexPrevPos(size_t i) const;

		size_t numNeiVertices(IdType iV) const;
		size_t neiVerticesStart(IdType iV) const;
		IdType getVertexIthNeiVertex(IdType iV, IdType neiId) const;

		size_t numNeiFaces(IdType iV) const;
		IdType getVertexIthNeiFace(IdType iV, IdType neiId) const;
		IdType getVertexIthNeiFaceOrder(IdType iV, IdType neiId) const;

		size_t numNeiEdges(IdType iV) const;
		IdType getVertexIthNeiEdge(IdType iV, IdType neiId) const;
		IdType getVertexIthNeiEdgeOrder(IdType iV, IdType neiId) const;

		size_t numRelevantBendings(IdType iV) const;
		IdType getVertexIthRelevantBending(IdType iV, IdType bendingId) const;
		IdType getVertexIthRelevantBendingOrder(IdType iV, IdType bendingId) const;
		

		void getEdgeVertexOrderInNeighborFace(int iE, int iNeighborFace, int& corner0, int& corner1) const {
			const EdgeInfo& edgeInfo = pTopology->edgeInfos[iE];
			if (0 != iNeighborFace && 1 != iNeighborFace) {
				assert(false);
			}

			corner0 = edgeInfo.eV1FaceOrder[iNeighborFace];
			corner1 = edgeInfo.eV2FaceOrder[iNeighborFace];
		}

		const EdgeInfo& getEdgeInfo(int eId) const
		{
			return pTopology->edgeInfos[eId];
		}

		std::vector<std::vector<int32_t>>& verticesColoringCategories() { return pTopology->verticesColoringCategories; }
		std::vector<int32_t> globalColors{};
		virtual int tearMesh(IdType v1, IdType v2);

	public:
		TriMeshParams::SharedPtr pObjectParams;


		// shape of positions is 3x(numVertices_ + 1), we pad an extra column for embree
		TVerticesMat positionsPrev;
		TVerticesMat velocities;
		TVerticesMat velocitiesPrev;
		TVerticesUVMat UVs;

		VecDynamic vertexMass;
		VecDynamic vertexInvMass;
		VecDynamic faceRestposeArea;

		VecDynamic DmInvs;

		std::vector<Eigen::Matrix<FloatingType, 2, 3>> restposeTriangles2D;

		FaceVIdsMat facePos;
		FaceVIdsMat faceUVs;
		TVerticesMat faceNormals;

		/*
		* constraints
		*/
		VecDynamicBool fixedMask;

		TriMeshTopology::SharedPtr pTopology;
		// records all the topology of the previously loaded mesh
		// we can grab reused the one that's already been computed since they don't change over time
		static std::map<std::string, TriMeshTopology::SharedPtr> topologies;
		static std::mutex topologies_lock;

		bool updated = true;
		bool activeForSim = true;
	private:
		int numVertices_;
		// do not access this directly, use positions() instead
		// because it's padded with an extra column for embree
		TVerticesMat positions_;
	};


	inline int TriMeshFEM::numFaces()
	{
		return facePos.cols();
	}

	inline int TriMeshFEM::numEdges()
	{
		assert(pTopology != nullptr);
		return pTopology->numEdges;
	}

	inline int TriMeshFEM::numVertices()
	{
		return numVertices_;
	}

	inline int TriMeshFEM::facePosVId(int tId, int vId) const 
	{
		return facePos(vId, tId);
	}

	inline Vec3Block TriMeshFEM::vertex(size_t i) 
	{
		return positions_.block<3, 1>(0, i);
	}

	inline ConstVec3Block TriMeshFEM::vertex(size_t i) const
	{
		return positions_.block<3, 1>(0, i);
	}

	inline ConstVec3Block TriMeshFEM::vertexPrevPos(size_t i) const
	{
		return positionsPrev.block<3, 1>(0, i);
	}

	inline size_t TriMeshFEM::neiVerticesStart(IdType iV) const
	{
		return pTopology->vertexNeighborVertices_infos(iV * 2);
	}

	inline IdType TriMeshFEM::getVertexIthNeiVertex(IdType iV, IdType neiId) const
	{
		return pTopology->vertexNeighborVertices(pTopology->vertexNeighborVertices_infos(iV * 2) + neiId);
	}

	inline size_t TriMeshFEM::numNeiFaces(IdType iV) const
	{
		return pTopology->vertexNeighborFaces_infos(iV * 2 + 1);
	}

	inline IdType TriMeshFEM::getVertexIthNeiFace(IdType iV, IdType neiId) const
	{
		return pTopology->vertexNeighborFaces(pTopology->vertexNeighborFaces_infos(iV * 2) + neiId);
	}

	inline IdType TriMeshFEM::getVertexIthNeiFaceOrder(IdType iV, IdType neiId) const
	{
		return pTopology->vertexNeighborFaces_vertexOrder(pTopology->vertexNeighborFaces_infos(iV * 2) + neiId);
	}

	inline size_t TriMeshFEM::numNeiEdges(IdType iV) const
	{
		return pTopology->vertexNeighborEdges_infos(iV*2+1);
	}

	inline IdType TriMeshFEM::getVertexIthNeiEdge(IdType iV, IdType neiId) const
	{
		return pTopology->vertexNeighborEdges(pTopology->vertexNeighborEdges_infos(iV * 2) + neiId);
	}

	inline IdType TriMeshFEM::getVertexIthNeiEdgeOrder(IdType iV, IdType neiId) const
	{
		return pTopology->vertexNeighborEdges_vertexOrder(pTopology->vertexNeighborEdges_infos(iV * 2) + neiId);
	}

	inline size_t TriMeshFEM::numRelevantBendings(IdType iV) const
	{
		return pTopology->vertexRelevantBendings_infos(iV*2+1);
	}

	inline IdType TriMeshFEM::getVertexIthRelevantBendingOrder(IdType iV, IdType bendingId) const
	{
		return pTopology->vertexRelevantBendings_vertexOrder(pTopology->vertexRelevantBendings_infos(iV * 2) + bendingId);
	}

	inline IdType TriMeshFEM::getVertexIthRelevantBending(IdType iV, IdType bendingId) const
	{
		return pTopology->vertexRelevantBendings(pTopology->vertexRelevantBendings_infos(iV * 2) + bendingId);
	}

	inline int TriMeshFEM::tearMesh(IdType v1, IdType v2)
	{
		int neiFaces1 = numNeiFaces(v1);
		int neiFaces2 = numNeiFaces(v2);
		int fId1 = -1, fId2 = -1; // fId1 aligns with (v1, v2), fId2 aligns with (v2, v1)
		for (int i = 0; i < neiFaces1; i++)
		{
			int fId = getVertexIthNeiFace(v1, i);
			if (facePosVId(fId, 0) == v1 && facePosVId(fId, 1) == v2 || facePosVId(fId, 1) == v1 && facePosVId(fId, 2) == v2 || facePosVId(fId, 2) == v1 && facePosVId(fId, 0) == v2)
			{
				fId1 = fId;
			}
			else if (facePosVId(fId, 0) == v2 && facePosVId(fId, 1) == v1 || facePosVId(fId, 1) == v2 && facePosVId(fId, 2) == v1 || facePosVId(fId, 2) == v2 && facePosVId(fId, 0) == v1)
			{
				fId2 = fId;
			}
		}
		if (fId1 == -1 || fId2 == -1)
		{
			std::cout << "Error: cannot find the two faces that share the edge (" << v1 << ", " << v2 << ")" << std::endl;
			return 0;
		}



		// change pTopology->faces3NeighborFaces of fId2
		if (facePosVId(fId2, 0) == v1)
		{
			pTopology->faces3NeighborFaces(2, fId2) = -1;
		}
		else if (facePosVId(fId2, 1) == v1)
		{
			pTopology->faces3NeighborFaces(0, fId2) = -1;
		}
		else if (facePosVId(fId2, 2) == v1)
		{
			pTopology->faces3NeighborFaces(1, fId2) = -1;
		}

		// change pTopology->faces3NeighborFaces of fId1
		if (facePosVId(fId1, 0) == v1) {
			pTopology->faces3NeighborFaces(0, fId1) = -1;
		}
		else if (facePosVId(fId1, 1) == v1) {
			pTopology->faces3NeighborFaces(1, fId1) = -1;
		}
		else if (facePosVId(fId1, 2) == v1) {
			pTopology->faces3NeighborFaces(2, fId1) = -1;
		}
		// v3 is the potential duplicate of v1
		// v4 is the potential duplicate of v2
		std::vector<int> faces1{ fId1 }; // faces belong to v1
		std::vector<int> faces2{ fId1 }; // faces belong to v2
		std::vector<int> faces3{ fId2 }; // faces belong to v3
		std::vector<int> faces4{ fId2 }; // faces belong to v4
		FloatingType area1 = faceRestposeArea(fId1);
		FloatingType area2 = faceRestposeArea(fId1);
		FloatingType area3 = faceRestposeArea(fId2);
		FloatingType area4 = faceRestposeArea(fId2);
		// rotate around v1 in two directions
		int fId = fId1;
		while (fId >= 0) {
			if (facePosVId(fId, 0) == v1) {
				fId = pTopology->faces3NeighborFaces(2, fId);
			}
			else if (facePosVId(fId, 1) == v1) {
				fId = pTopology->faces3NeighborFaces(0, fId);
			}
			else if (facePosVId(fId, 2) == v1) {
				fId = pTopology->faces3NeighborFaces(1, fId);
			}
			if (fId >= 0) {
				faces1.push_back(fId);
				area1 += faceRestposeArea(fId);
			}
		}
		fId = fId2;
		while (fId >= 0) {
			if (facePosVId(fId, 0) == v1) {
				fId = pTopology->faces3NeighborFaces(0, fId);
			}
			else if (facePosVId(fId, 1) == v1) {
				fId = pTopology->faces3NeighborFaces(1, fId);
			}
			else if (facePosVId(fId, 2) == v1) {
				fId = pTopology->faces3NeighborFaces(2, fId);
			}
			if (fId >= 0) {
				faces3.push_back(fId);
				area3 += faceRestposeArea(fId);
			}
		}
		// rotate around v2 in two directions
		fId = fId1;
		while (fId >= 0) {
			if (facePosVId(fId, 0) == v2) {
				fId = pTopology->faces3NeighborFaces(0, fId);
			}
			else if (facePosVId(fId, 1) == v2) {
				fId = pTopology->faces3NeighborFaces(1, fId);
			}
			else if (facePosVId(fId, 2) == v2) {
				fId = pTopology->faces3NeighborFaces(2, fId);
			}
			if (fId >= 0) {
				faces2.push_back(fId);
				area2 += faceRestposeArea(fId);
			}
		}
		fId = fId2;
		while (fId >= 0) {
			if (facePosVId(fId, 0) == v2) {
				fId = pTopology->faces3NeighborFaces(2, fId);
			}
			else if (facePosVId(fId, 1) == v2) {
				fId = pTopology->faces3NeighborFaces(0, fId);
			}
			else if (facePosVId(fId, 2) == v2) {
				fId = pTopology->faces3NeighborFaces(1, fId);
			}
			if (fId >= 0) {
				faces4.push_back(fId);
				area4 += faceRestposeArea(fId);
			}
		}
		int ret = 0;
		if (faces1.size() + faces3.size() == neiFaces1) {
			// new vertex
			ret |= 1;
			int v3 = numVertices_;
			numVertices_++;
			positions_.conservativeResize(Eigen::NoChange, numVertices_);
			positions_.col(v3) = positions_.col(v1);
			positionsPrev.conservativeResize(Eigen::NoChange, numVertices_);
			positionsPrev.col(v3) = positionsPrev.col(v1);
			velocities.conservativeResize(Eigen::NoChange, numVertices_);
			velocities.col(v3) = velocities.col(v1);
			velocitiesPrev.conservativeResize(Eigen::NoChange, numVertices_);
			velocitiesPrev.col(v3) = velocitiesPrev.col(v1);
			fixedMask.conservativeResize(numVertices_);
			fixedMask(v3) = fixedMask(v1);

			// change mass later
			vertexMass.conservativeResize(numVertices_);
			vertexInvMass.conservativeResize(numVertices_);

			// update facePos
			for (int i = 0; i < faces3.size(); ++i) {
				if (facePosVId(faces3[i], 0) == v1) {
					facePos(0, faces3[i]) = v3;
				}
				else if (facePosVId(faces3[i], 1) == v1) {
					facePos(1, faces3[i]) = v3;
				}
				else if (facePosVId(faces3[i], 2) == v1) {
					facePos(2, faces3[i]) = v3;
				}
			}

			// update pTopology->vertexNeighborFaces and pTopology->vertexNeighborFaces_infos
			pTopology->vertexNeighborFaces_infos.conservativeResize(numVertices_ * 2);
			pTopology->vertexNeighborFaces_infos(v1 * 2 + 1) = faces1.size();
			pTopology->vertexNeighborFaces_infos(v3 * 2 + 1) = faces3.size();
			pTopology->vertexNeighborFaces_infos(v3 * 2) = pTopology->vertexNeighborFaces_infos(v1 * 2) + faces1.size();
			for (int i = 0; i < faces1.size(); ++i) {
				pTopology->vertexNeighborFaces(pTopology->vertexNeighborFaces_infos(v1 * 2) + i) = faces1[i];
			}
			for (int i = 0; i < faces3.size(); ++i) {
				pTopology->vertexNeighborFaces(pTopology->vertexNeighborFaces_infos(v3 * 2) + i) = faces3[i];
			}

			// update mass
			vertexMass(v3) = vertexMass(v1) * area3 / (area1 + area3);
			vertexMass(v1) = vertexMass(v1) * area1 / (area1 + area3);
			vertexInvMass(v3) = 1 / vertexMass(v3);
			vertexInvMass(v1) = 1 / vertexMass(v1);

			// inherit color
			globalColors.resize(numVertices_);
			globalColors[v3] = globalColors[v1];
		}
		else {
			assert(faces1.size() == neiFaces1 && faces3.size() == neiFaces1);
		}
		if (faces2.size() + faces4.size() == neiFaces2) {
			// new vertex
			ret |= 2;
			int v4 = numVertices_;
			numVertices_++;
			positions_.conservativeResize(Eigen::NoChange, numVertices_);
			positions_.col(v4) = positions_.col(v2);
			positionsPrev.conservativeResize(Eigen::NoChange, numVertices_);
			positionsPrev.col(v4) = positionsPrev.col(v2);
			velocities.conservativeResize(Eigen::NoChange, numVertices_);
			velocities.col(v4) = velocities.col(v2);
			velocitiesPrev.conservativeResize(Eigen::NoChange, numVertices_);
			velocitiesPrev.col(v4) = velocitiesPrev.col(v2);
			fixedMask.conservativeResize(numVertices_);
			fixedMask(v4) = fixedMask(v2);

			// change mass later
			vertexMass.conservativeResize(numVertices_);
			vertexInvMass.conservativeResize(numVertices_);

			// update facePos
			for (int i = 0; i < faces4.size(); ++i) {
				if (facePosVId(faces4[i], 0) == v2) {
					facePos(0, faces4[i]) = v4;
				}
				else if (facePosVId(faces4[i], 1) == v2) {
					facePos(1, faces4[i]) = v4;
				}
				else if (facePosVId(faces4[i], 2) == v2) {
					facePos(2, faces4[i]) = v4;
				}
			}

			// update pTopology->vertexNeighborFaces and pTopology->vertexNeighborFaces_infos
			pTopology->vertexNeighborFaces_infos.conservativeResize(numVertices_ * 2);
			pTopology->vertexNeighborFaces_infos(v2 * 2 + 1) = faces2.size();
			pTopology->vertexNeighborFaces_infos(v4 * 2 + 1) = faces4.size();
			pTopology->vertexNeighborFaces_infos(v4 * 2) = pTopology->vertexNeighborFaces_infos(v2 * 2) + faces2.size();
			for (int i = 0; i < faces2.size(); ++i) {
				pTopology->vertexNeighborFaces(pTopology->vertexNeighborFaces_infos(v2 * 2) + i) = faces2[i];
			}
			for (int i = 0; i < faces4.size(); ++i) {
				pTopology->vertexNeighborFaces(pTopology->vertexNeighborFaces_infos(v4 * 2) + i) = faces4[i];
			}

			// update mass
			vertexMass(v4) = vertexMass(v2) * area4 / (area2 + area4);
			vertexMass(v2) = vertexMass(v2) * area2 / (area2 + area4);
			vertexInvMass(v4) = 1 / vertexMass(v4);
			vertexInvMass(v2) = 1 / vertexMass(v2);

			// inherit color
			globalColors.resize(numVertices_);
			globalColors[v4] = globalColors[v2];
		}
		else {
			assert(faces2.size() == neiFaces2 && faces4.size() == neiFaces2);
		}
		return ret;
	}

	inline size_t TriMeshFEM::numNeiVertices(IdType iV) const 
	{
		return pTopology->vertexNeighborVertices_infos(iV * 2 + 1);
	}
	 
	inline Vec3 TriMeshFEM::computeNormal(int fId, bool normalize) const
	{
		Vec3 v1 = positions_.col(facePosVId(fId, 1)) - positions_.col(facePosVId(fId, 0));
		Vec3 v2 = positions_.col(facePosVId(fId, 2)) - positions_.col(facePosVId(fId, 0));

		Vec3 normal = v1.cross(v2);
		if (normalize)
		{
			normal.normalize();

		}

		return normal;
	}

	inline Eigen::Map<Mat2> TriMeshFEM::getDmInv(int fId)
	{
		return Eigen::Map<Mat2>(DmInvs.data() + fId * 4);
	}

	inline Eigen::Block<TVerticesMat, 3, -1> TriMeshFEM::positions()
	{
		return positions_.block<3, -1>(0, 0, 3, numVertices_);
	}
}