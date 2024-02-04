#include "TetMeshFEM.h"
#include "../Parallelization/CPUParallelization.h"

#include <MeshFrame/Utility/IO.h>
#include <random>
using namespace GAIA;
#include <MeshFrame/Utility/Parser.h>
#include <chrono>       // std::chrono::system_clock
#include <MeshFrame/Memory/Array.h>
#include <unordered_set>

#include "../IO/FileIO.h"

#define SIZE_CANDIDATE_FACE_STACK 32
#define SIZE_TRAVERSED_LIST_STACK 128
#define SIZE_TRAVERSED_CIRCULAR_ARRAY 32
// #define TRAVERESED_TETS_SET
#define TRAVERESED_TETS_CIRCULAR_ARRAY
//#define ONLY_RECORD_BRANCHING_TETS

const int32_t TetMeshFEM::passThroughCheckSteps = 5;
const int32_t TetMeshFEM::tet4Faces[4][3] = { { 1, 2, 3 },{ 2, 0, 3 },{ 0, 1, 3 },{ 1, 0, 2 } };
// how the 3 possible exit face is made of; note that the indices corresponds to ptsProj2D 
// (which have been re-permuted such that the incoming face will be f3
const int32_t TetMeshFEM::posibleExit3Faces[3][3] = { { 1, 2, 3 },{ 0, 2, 3 },{ 0, 1, 3 } };
std::map<std::string, TetMeshTopology::SharedPtr> TetMeshFEM::topologies;
std::mutex TetMeshFEM::topologies_lock;

//Vec3& GAIA::TetMesh::getV(int vId)
//{
//	auto& v = vertices.block<3, 1>(0, vId);
//
//	return v;
namespace GAIA {
	typedef MF::TriMesh::CIterators<TetSurfaceMeshMF> ItSurface;

}


void GAIA::TetMeshTopology::initialize(TetMeshMF * pTM_MF, ObjectParams::SharedPtr pObjectParams)
{
	nVerts = pTM_MF->numVertices();
	tetVIds = pTM_MF->tetVIds();
	// initialize topological data
	// - tet's vertex XOR sum and neighbor tets
	tetsXorSums.resize(numTets());
	tetsNeighborTets.resizeLike(tetVIds);
	for (int iTet = 0; iTet < numTets(); iTet++)
	{
		tetsXorSums(iTet) = tetVIds(0, iTet) ^ tetVIds(1, iTet) ^ tetVIds(2, iTet) ^ tetVIds(3, iTet);

		// const int32_t order[4][3] = { { 1, 2, 3 },{ 2, 0, 3 },{ 0, 1, 3 },{ 1, 0, 2 } };
		for (int iV = 0; iV < 4; iV++)
		{
			/*int32_t vId = tetsXorSums(iTet) ^ tetVIds(order[iV][0], iTet) ^ tetVIds(order[iV][1], iTet) ^ tetVIds(order[iV][2], iTet);
			assert(vId == tetVIds(iV, iTet));*/

			// - compute tet's 4 neighbor tets
			assert(tetVIds(iV, iTet) == pTM_MF->TetTVertex(pTM_MF->idTet(iTet), iV)->vert()->id());
			TetMeshMF::HFPtr pOppositeHF = pTM_MF->TVertexOppositeHalfFace(pTM_MF->TetTVertex(pTM_MF->idTet(iTet), iV))->dual();

			int32_t neiTetId = -1;
			if (pOppositeHF != nullptr)
			{
				neiTetId = pOppositeHF->tet()->id();
			}

			tetsNeighborTets(iV, iTet) = neiTetId;
		}
	}

#ifdef TET_TET_ADJACENT_LIST
	tetAllNeighborTets.resize(tetVIds.cols());
	cpu_parallel_for(0, numTets(), [&](int iTet) {
		std::unordered_set<int> neighborTetsId;
		TetMeshMF::TPtr pT = pTM_MF->idTet(iTet);
		for (TetMeshMF::VPtr pV : TIt::T_VIterator(pT))
		{
			for (TetMeshMF::TVPtr pNeiTV : TIt::V_TVIterator(pV))
			{
				TetMeshMF::TPtr pNeiT = TetMeshMF::TVertexTet(pNeiTV);
				if (pNeiT != pT) {
					neighborTetsId.insert(pNeiT->id());
				}
			}
		}

		for (int id : neighborTetsId)
		{
			tetAllNeighborTets[iTet].push_back(id);
		}
		});
#endif // TET_TET_ADJACENT_LIST


	// load tet coloring information
	if (pObjectParams->tetsColoringCategoriesPath != "")
	{
		nlohmann::json tetsColoring;

		MF::loadJson(pObjectParams->tetsColoringCategoriesPath, tetsColoring);
		MF::convertJsonParameters(tetsColoring, tetsColoringCategories);

		// we sort them from large to small for the aggregated solve
		std::sort(tetsColoringCategories.begin(), tetsColoringCategories.end(),
			[](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });
	}

	// load edge coloring information
	if (pObjectParams->edgesColoringCategoriesPath != "")
	{
		nlohmann::json edgesColoring;

		MF::loadJson(pObjectParams->edgesColoringCategoriesPath, edgesColoring);
		MF::convertJsonParameters(edgesColoring, edgesColoringCategories);

		// we sort them from large to small for the aggregated solve
		std::sort(edgesColoringCategories.begin(), edgesColoringCategories.end(),
			[](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });
	}

	// load vertices coloring information
	if (pObjectParams->verticesColoringCategoriesPath != "")
	{
		nlohmann::json vertsColoring;

		MF::loadJson(pObjectParams->verticesColoringCategoriesPath, vertsColoring);
		MF::convertJsonParameters(vertsColoring, verticesColoringCategories);

		// we sort them from large to small for the aggregated solve
		std::sort(verticesColoringCategories.begin(), verticesColoringCategories.end(),
			[](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });
	}


	TetSurfaceMeshMF* pSurfaceMesh = pTM_MF->createSurfaceMesh<TetSurfaceMeshMF>();
	pSurfaceMesh->reinitializeVId();
	pSurfaceMesh->reinitializeFId();
	surfaceVIds.resize(pSurfaceMesh->numVertices());
	tetVertIndicesToSurfaceVertIndices.resize(pTM_MF->numVertices());
	tetVertIndicesToSurfaceVertIndices = decltype(tetVertIndicesToSurfaceVertIndices)::Constant(tetVertIndicesToSurfaceVertIndices.rows(),
		tetVertIndicesToSurfaceVertIndices.cols(), -1);
	surfaceVertexNeighborSurfaceFaces.resize(pSurfaceMesh->numVertices());
	surfaceVertexNeighborSurfaceVertices.resize(pSurfaceMesh->numVertices());
	tetsIsSurfaceTet = VecDynamicBool::Constant(numTets(), true);

	int iV = 0;
	for (TetSurfaceMeshMF::VPtr pSurfV : ItSurface::MVIterator(pSurfaceMesh))
	{
		assert(pSurfV->id() == iV);
		surfaceVIds(iV) = pSurfV->getTetMeshVertId();
		tetVertIndicesToSurfaceVertIndices(pSurfV->getTetMeshVertId()) = iV;

		for (TetSurfaceMeshMF::VPtr pSurfVNei : ItSurface::VVIterator(pSurfV)) {
			surfaceVertexNeighborSurfaceVertices[iV].push_back(pSurfVNei->getTetMeshVertId());
		}

		for (TetSurfaceMeshMF::FPtr pSurfFNei : ItSurface::VFIterator(pSurfV)) {
			surfaceVertexNeighborSurfaceFaces[iV].push_back(pSurfFNei->id());
		}

		TetMeshMF::VPtr pV = pTM_MF->idVertex(pSurfV->getTetMeshVertId());
		for (TetMeshMF::TVPtr pNeiTV : TIt::V_TVIterator(pV))
		{
			TetMeshMF::TPtr pNeiT = TetMeshMF::TVertexTet(pNeiTV);
			tetsIsSurfaceTet[pNeiT->id()] = true;
		}
		iV++;
	}

	surfaceFacesTetMeshVIds.resize(3, pSurfaceMesh->numFaces());
	surfaceFacesSurfaceMeshVIds.resize(3, pSurfaceMesh->numFaces());
	surfaceFacesBelongingTets.resize(pSurfaceMesh->numFaces());
	surfaceFacesIdAtBelongingTets.resize(pSurfaceMesh->numFaces());
	surfaceFaces3NeighborFaces.resize(3, pSurfaceMesh->numFaces());

	int iF = 0;
	for (TetSurfaceMeshMF::FPtr pSurfF : ItSurface::MFIterator(pSurfaceMesh))
	{
		assert(pSurfF->id() == iF);

		size_t iV = 0;
		for (TetSurfaceMeshMF::HEPtr pHE : ItSurface::FHEIterator(pSurfF))
		{
			surfaceFacesTetMeshVIds(iV, iF) = TetSurfaceMeshMF::halfedgeSource(pHE)->getTetMeshVertId();
			// say 3 vertices in surfaceFacesTetMeshVIds are A, B and C,
			// the 3 edges will be AB, BC, CD
			// and 3 neighbor faces will be on three face on the other side of AB, BC, CA correspondingly
			// this gurrantees that
			surfaceFaces3NeighborFaces(iV, iF) = TetSurfaceMeshMF::halfedgeSym(pHE)->face()->id();
			++iV;
		}

		surfaceFacesBelongingTets(iF) = pSurfF->getTetMeshHalfFacePtr()->tet()->id();
		// face id is the vId that is not included by this face
		surfaceFacesIdAtBelongingTets(iF) = -1;
		for (int iFTet = 0; iFTet < 4; iFTet++)
		{
			int32_t vId = getTVId(surfaceFacesBelongingTets(iF), iFTet);
			bool faceHasVId = false;
			for (int fvId = 0; fvId < 3; fvId++)
			{
				if (surfaceFacesTetMeshVIds(fvId, iF) == vId)
				{
					faceHasVId = true;
					break;
				}
			}
			if (!faceHasVId)
			{
				surfaceFacesIdAtBelongingTets(iF) = iFTet;
			}
		}
		assert(surfaceFacesIdAtBelongingTets(iF) != -1);

		surfaceFacesSurfaceMeshVIds(0, iF) = tetVertIndicesToSurfaceVertIndices(surfaceFacesTetMeshVIds(0, iF));
		surfaceFacesSurfaceMeshVIds(1, iF) = tetVertIndicesToSurfaceVertIndices(surfaceFacesTetMeshVIds(1, iF));
		surfaceFacesSurfaceMeshVIds(2, iF) = tetVertIndicesToSurfaceVertIndices(surfaceFacesTetMeshVIds(2, iF));

		++iF;
	}

	// intialize Edges
	nEdges = pTM_MF->numEdges();
	edges.resize(2, nEdges);

	int iEdge = 0;
	std::vector<IdType> edgeNeighborTets_;
	std::vector<IdType> edgeNighborTets_Vert1TetOrders_;
	std::vector<IdType> edgeNighborTets_Vert2TetOrders_;
	edgeNeighborTets_infos.resize(2 * nEdges);

	for (TetMeshMF::EPtr pE : TIt::TM_EIterator(pTM_MF))
	{
		edges.col(iEdge) << pE->vertex1()->id(), pE->vertex2()->id();

		int32_t* edge = edges.data() + 2 * iEdge;
		assert(pE->vertex1()->id() == edge[0] && pE->vertex2()->id() == edge[1]);


		// compute neighbor tets
		size_t numNeiTets = 0;
		edgeNeighborTets_infos(iEdge * 2) = edgeNeighborTets_.size();
		for (TetMeshMF::TEPtr pENeiTE : TIt::E_TEIterator(pE))
		{
			TetMeshMF::TPtr pENeiT = TetMeshMF::TEdgeTet(pENeiTE);
			IdType tId = pENeiT->id();
			edgeNeighborTets_.push_back(tId);
			numNeiTets++;

			int eVId1 = pE->vertex1()->id();
			int eVId2 = pE->vertex2()->id();
			for (size_t iTVId = 0; iTVId < 4; iTVId++)
			{
				if (tetVIds(iTVId, tId) == eVId1) {
					edgeNighborTets_Vert1TetOrders_.push_back(iTVId);
				}

				if (tetVIds(iTVId, tId) == eVId2) {
					edgeNighborTets_Vert2TetOrders_.push_back(iTVId);
				}
			}
			assert(edgeNighborTets_Vert1TetOrders_.size() == edgeNighborTets_Vert2TetOrders_.size());

		}
		edgeNeighborTets_infos(iEdge * 2 + 1) = numNeiTets;
		++iEdge;
	}
	//std::cout << "edgeNeighborTets_infos:\n" << edgeNeighborTets_infos.transpose();
	edgeNeighborTets = VecDynamicI::Map(&edgeNeighborTets_[0], edgeNeighborTets_.size());
	edgeNeighborTets_edge2VerticesOrderInTet.resize(2, edgeNeighborTets.size());
	edgeNeighborTets_edge2VerticesOrderInTet.row(0) = VecDynamicI::Map(&edgeNighborTets_Vert1TetOrders_[0], edgeNighborTets_Vert1TetOrders_.size());
	edgeNeighborTets_edge2VerticesOrderInTet.row(1) = VecDynamicI::Map(&edgeNighborTets_Vert2TetOrders_[0], edgeNighborTets_Vert2TetOrders_.size());

	std::vector<IdType> vertexNeighborTets_;
	std::vector<IdType> vertexNeighborTets_tetVId_;

	iV = 0;

	size_t maxNumberOfNeiTets = 0;
	float avgNumberOfNeiTets = 0;

	vertexNeighborTets_infos.resize(2 * pTM_MF->numVertices());
	for (TetMeshMF::VPtr pV : TIt::TM_VIterator(pTM_MF))
	{
		size_t numNeiTets = 0;
		vertexNeighborTets_infos(iV * 2) = vertexNeighborTets_.size();
		assert(iV == pV->id());

		for (TetMeshMF::TVPtr pTV : TIt::V_TVIterator(pV))
		{
			++numNeiTets;
			TetMeshMF::TPtr pT = (TetMeshMF::TPtr)pTV->tet();
			IdType tId = pT->id();
			vertexNeighborTets_.push_back(tId);

			// compute tet VId

			int tetVId = -1;
			for (size_t iTVId = 0; iTVId < 4; iTVId++)
			{
				if (tetVIds(iTVId, tId) == iV) {
					tetVId = iTVId;
					break;
				}
			}
			assert(tetVId != -1);

			vertexNeighborTets_tetVId_.push_back(tetVId);
		}
		vertexNeighborTets_infos(iV * 2 + 1) = numNeiTets;

		avgNumberOfNeiTets += numNeiTets;
		if (maxNumberOfNeiTets < numNeiTets)
		{
			maxNumberOfNeiTets = numNeiTets;
		}

		iV++;
	}

	avgNumberOfNeiTets /= pTM_MF->numVertices();

	vertexNeighborTets = VecDynamicI::Map(&vertexNeighborTets_[0], vertexNeighborTets_.size());
	vertexNeighborTets_vertexOrder = VecDynamicI::Map(&vertexNeighborTets_tetVId_[0], vertexNeighborTets_tetVId_.size());

	// std::cout << "vertexNeighborTets: " << vertexNeighborTets.transpose() << "\n";
	// std::cout << "vertexNeighborTets_vertexOrder: " << vertexNeighborTets_vertexOrder.transpose() << "\n";
	// std::cout << "vertexNeighborTets_infos: " << vertexNeighborTets_infos.transpose() << "\n";
}

void GAIA::TetMeshFEM::initialize(ObjectParams::SharedPtr inObjectParams, std::shared_ptr<TetMeshMF> pTM_MF)
{
#ifdef KEEP_MESHFRAME_MESHES
	m_pTM_MF = pTM_MF;
#endif // KEEP_MESHFRAME_MESHES

	MF::IO::FileParts fp = MF::IO::fileparts(inObjectParams->path);

	//bool loadSucceed = false;
	//if (fp.ext == ".t")
	//{
	//	//pTM->_load_t(param.inputModelPath.c_str(), true);
	//	load_t(inObjectParams->path.c_str());
	//	loadSucceed = true;
	//}
	//else if (fp.ext == ".vtk") {
	//	std::cout << "Currently vtk file is not supported! " << std::endl;
	//	//continue;
	//}
	//else
	//{
	//	std::cout << "Unsupported file format: " << fp.ext << std::endl;
	//	//continue;
	//}
	// initialize the vertices and tets
	// we add one more vertex at the end to make sure it can be used as the vertex buffer of embree
	mVertPos = TVerticesMat::Zero(POINT_VEC_DIMS, pTM_MF->numVertices() + 1);
	mVertPos.block(0, 0, POINT_VEC_DIMS, pTM_MF->numVertices()) = pTM_MF->vertPos();
	//TBB_PARALLEL_FOR(0, tetVIds.cols(), iTet, );

	m_nVertices = pTM_MF->numVertices();
	m_nTets = pTM_MF->numTets();
	pObjectParams = inObjectParams;

	computeTopology(pTM_MF.get());


	Eigen::AngleAxisf angleAxis;
	FloatingType rot = pObjectParams->rotation.norm();
	angleAxis.angle() = rot;
	angleAxis.axis() = pObjectParams->rotation / rot;

	if (rot > 0)
	{
		Eigen::Matrix3f R;
		R = angleAxis.toRotationMatrix();
		mVertPos = R * mVertPos;
	}

	mVertPos.row(0) *= pObjectParams->scale(0);
	mVertPos.row(1) *= pObjectParams->scale(1);
	mVertPos.row(2) *= pObjectParams->scale(2);
	mVertPos.colwise() += pObjectParams->translation;

	vertexMass = decltype(vertexMass)::Zero(m_nVertices);

	tetRestVolume.resize(m_nTets);
	tetInvRestVolume.resize(m_nTets);
	DmInvs.resize(m_nTets * 9);

	verticesInvertedSign.resize(m_nVertices);
	verticesInvertedSign.setZero();
	verticesInvertedSignPrevPos.resize(m_nVertices);
	verticesInvertedSignPrevPos.setZero();

	tetsInvertedSign.resize(m_nTets);
	tetsInvertedSign.setZero();
	tetsInvertedSignPrevPos.resize(m_nTets);
	tetsInvertedSignPrevPos.setZero();

	verticesCollisionDetectionEnabled = VecDynamicBool::Constant(numVertices(), true);

	mVelocity = TVerticesMat::Zero(mVertPos.rows(), mVertPos.cols());
	mVelocity.colwise() = pObjectParams->initialVelocity;

	mVertPrevPos = TVerticesMat::Zero(mVertPos.rows(), mVertPos.cols());

	for (int iTet = 0; iTet < numTets(); ++iTet)
	{
		Mat3 Dm;
		Eigen::Map<Mat3> DmInv = getDmInv(iTet);
		computeDs(Dm, iTet);
		//std::cout << "Ds:\n" << Ds << "\n";
		DmInv = Dm.inverse();
		//std::cout << "DsInv:\n" << DSInv << "\n";

		// compute volume and mass
		FloatingType vol = Dm.determinant() / 6;
		// compute DsInvs;

		tetRestVolume[iTet] = vol;
		tetInvRestVolume[iTet] = 1. / vol;
		for (size_t iV = 0; iV < 4; iV++)
		{
			vertexMass(pTopology->tetVIds(iV, iTet)) += 0.25 * vol * pObjectParams->density;
		}
	}
	vertexInvMass = vertexMass.cwiseInverse();

#ifdef ENABLE_REST_POSE_CLOSEST_POINT
	restposeVerts = mVertPos;
#endif // ENABLE_REST_POSE_CLOSEST_POINT

	fixedMask.resize(numVertices());
	fixedMask.setZero();

	for (size_t iV = 0; iV < pObjectParams->fixedPoints.size(); iV++)
	{
		fixedMask(pObjectParams->fixedPoints[iV]) = true;
	}

}

size_t GAIA::TetMeshFEM::numVertices()
{
	return m_nVertices;
}

size_t GAIA::TetMeshFEM::numTets()
{
	return m_nTets;
}

void GAIA::TetMeshFEM::save(const std::string outFile)
{
	MF::IO::FileParts fp(outFile);

	if (fp.ext == ".ply") {
		saveAsPLY(outFile);
	}
	else if (fp.ext == ".obj")
	{
		std::cout << "Warning! saveToOBJ has not be implemented yet!\n";
	}
	else
	{
		std::cout << "Error!!! Unrecognized extension name:" << fp.ext << "!\n";
	}

}

void GAIA::TetMeshFEM::saveAsPLY(const std::string outFile)
{
	std::vector<std::array<FloatingType, 3>> verts;
	std::vector<std::array<int, 3>> faces;

	for (int iV = 0; iV < surfaceVIds().size(); ++iV)
	{
		auto v = surfaceVertex(iV);
		std::array<FloatingType, 3> pt = { v[0], v[1], v[2] };

		verts.push_back(pt);
	}

	for (int iF = 0; iF < surfaceFacesSurfaceMeshVIds().cols(); ++iF)
	{
		std::array<int, 3> fVIds;
		for (size_t iFV = 0; iFV < 3; iFV++)
		{
			fVIds[iFV] = surfaceFacesSurfaceMeshVIds()(iFV, iF);
		}
		faces.push_back(fVIds);
	}

	saveAsPly(outFile, verts, faces);

	return;
}



bool TetMeshFEM::tetrahedralTraverseTo(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
	int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics)
{
	statistics.stopReason = TraverseStopReason::querySuccess;

	if (startTetId == targetTetId)
	{
		return true;
	}
	// get the ray-aligned 2D coordinate system
	Eigen::Matrix<FloatingType, 3, 2> axes;
	CuMatrix::buildOrthonormalBasis(rayDir.data(), axes.col(0).data(), axes.col(1).data());
	
	assert(axes.col(0).dot(axes.col(1)) < 1e-6);
	assert(axes.col(0).dot(rayDir) < 1e-6);
	assert(axes.col(1).dot(rayDir) < 1e-6);

	Eigen::Matrix<FloatingType, 2, 3> axesT = axes.transpose();

	// initialize the 2D coordinate

	Eigen::Matrix<FloatingType, 2, 4> ptsProj2D;

	int count = 0;
	// rearrange the ording of vertices thus the vertex across the incoming face is the last one
	// and the first 3 vertices are ordered as tet4Faces[startFaceId]
	projectTo2DCoordinates(ptsProj2D, startTetId, startFaceId, axesT, rayOrigin);

	Eigen::Vector4i possibleExitFace;
	exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);

	// figure out the outcoming Face
	CPArrayStatic<int32_t, SIZE_CANDIDATE_FACE_STACK> candidateExitFaces;
	CPArrayStatic<int32_t, SIZE_CANDIDATE_FACE_STACK> candidateExitFacesCorrespondingTets;

#if defined TRAVERESED_TETS_SET
	std::unordered_set<int32_t> traversedTets;
	traversedTets.insert(startTetId);
#elif defined TRAVERESED_TETS_CIRCULAR_ARRAY
	CircularArray<int32_t, SIZE_TRAVERSED_CIRCULAR_ARRAY> traversedTets;
	traversedTets.push_back(startTetId);
#else
	CPArrayStatic<int32_t, SIZE_TRAVERSED_LIST_STACK> traversedTets;
	traversedTets.push_back(startTetId);
#endif

#ifdef OUTPUT_TRAVERSED_TETS 
	statistics.traversedTetsOutput.clear();
	statistics.traversedTetsOutput.push_back(startTetId);
#endif

	for (int i = 0; i < 3; i++)
	{
		if (possibleExitFace(i))
		{
			// the first 3 vertices are ordered as tet4Faces[startFaceId]
			// thus tet4Faces[startFaceId][i] is the exit face id
			int32_t exitFaceId = tet4Faces[startFaceId][i];

			candidateExitFaces.push_back(exitFaceId);
			candidateExitFacesCorrespondingTets.push_back(startTetId);
		}
	}

	statistics.numTetsTraversed = 0;
	while (!candidateExitFaces.empty())
	{
		int32_t intersectedFaceId = candidateExitFaces.pop_back();
		int32_t previousTetId = candidateExitFacesCorrespondingTets.pop_back();

		int32_t intersectedFaceVIds[3] = {
			tetVIds()(tet4Faces[intersectedFaceId][0], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][1], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][2], previousTetId)
		};

		int32_t currentTetId = tetsNeighborTets()(intersectedFaceId, previousTetId);

		//std::cout << "Previous tet VIds: " << tetVIds.col(previousTetId).transpose() << std::endl;

#ifdef OUTPUT_TRAVERSED_TETS 
		statistics.traversedTetsOutput.push_back(currentTetId);
#endif

		if (currentTetId == targetTetId)
		{
			statistics.stopReason = TraverseStopReason::querySuccess;
			return true;
		}

		if (currentTetId == -1)
			// reached the boundary
		{
			statistics.stopReason = TraverseStopReason::reachedBoundary;
			return false;
		}
#if defined TRAVERESED_TETS_SET

		if (traversedTets.find(currentTetId) != traversedTets.end())
			// formed loop, cut this branch here
		{
			continue;
		}
		else
		{
			traversedTets.insert(currentTetId);
		}
#elif defined TRAVERESED_TETS_CIRCULAR_ARRAY
		if (traversedTets.hasBackward(currentTetId))
		{
			continue;
		}
#ifndef ONLY_RECORD_BRANCHING_TETS
		traversedTets.push_back(currentTetId);
#endif // ONLY_RECORD_BRANCHING_TETS

#else
		if (traversedTets.has(currentTetId))
			// formed loop, cut this branch here
		{
			continue;
		}
		else
			// add currentTetId to traversed tets list
		{
			if (!traversedTets.push_back(currentTetId)) {
				stopReason = TraverseStopReason::overflow;
				return false;
			}
		}
#endif

		int32_t newVertexId = tetsXorSums()(currentTetId) 
			^ intersectedFaceVIds[0]
			^ intersectedFaceVIds[1]
			^ intersectedFaceVIds[2];
		
		// the intersectedFaceId is not the incoming face id of intersectedFaceId
		// we need to query the id of the incoming face
		int32_t incomingFaceIdCurTet = 3;
		for (int iV = 0; iV < 3; iV++)
		{
			if (newVertexId == tetVIds()(iV, currentTetId)) {
				incomingFaceIdCurTet = iV;

				break;
			}
		}
		//std::cout << "Current tet VIds: " << tetVIds.col(currentTetId).transpose() << std::endl;

		// rearrange the ording of vertices thus the vertex across the incoming face (incomingFaceIdCurTet) is the last one
		// and the first 3 vertices are ordered as tet4Faces[startFaceId]
		projectTo2DCoordinates(ptsProj2D, currentTetId, incomingFaceIdCurTet, axesT, rayOrigin);
		int numNumIntersectingFaces = exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);
#ifdef ONLY_RECORD_BRANCHING_TETS
		if (numNumIntersectingFaces > 1)
		{
			traversedTets.push_back(currentTetId);
		}
#endif // ONLY_RECORD_BRANCHING_TETS

		for (int iF = 0; iF < 3; iF++)
		{
			if (possibleExitFace(iF))
				// add intersected faces to candidate intersection face stack
			{
				// the exit face id is the actual vertex id (0~3), not how they are arranged in ptsProj2D
				// exit face id can be deduced from tet4Faces[incomingFaceIdCurTet]
				int32_t exitFaceId = tet4Faces[incomingFaceIdCurTet][iF];
				int32_t nextTetId = tetsNeighborTets()(exitFaceId, currentTetId);
				if (nextTetId == targetTetId)
				{
					statistics.stopReason = TraverseStopReason::querySuccess;
					return true;
				}

				if (!candidateExitFaces.push_back(exitFaceId) || !candidateExitFacesCorrespondingTets.push_back(currentTetId))
					// the static array has overflown, static version of this function cannot be used
					// in this case this function will return false and set stackOverflow=true, to notify the frame work to call a version with dynamic mem allocation
				{
					statistics.stopReason = TraverseStopReason::overflow;
					return false;
				}

				if ((!(statistics.numTetsTraversed % passThroughCheckSteps)) && maxTraversalDis > 0.f)
					// check pass through
				{
					Vec3 barys, intersectPt;


					//std::cout << "ptsProj2D: \n" << ptsProj2D << std::endl;

					originBarycentric2DTriangle(ptsProj2D.col(posibleExit3Faces[iF][0]),
						ptsProj2D.col(posibleExit3Faces[iF][1]),
						ptsProj2D.col(posibleExit3Faces[iF][2]),
						barys);
					intersectPt = Vec3::Zero();

					//for (int iV = 0; iV < 3; iV++)
					//{
					//	intersectPt += barys(iV) * ptsPermuted3D.col(posibleExit3Faces[iF][iV]);
					//}
					Eigen::Matrix<FloatingType, 3, 4> ptsPermuted3D;
					copyRepermutedVerts(ptsPermuted3D, currentTetId, incomingFaceIdCurTet);

					intersectPt += barys(0) * ptsPermuted3D.col(posibleExit3Faces[iF][0])
					            +  barys(1) * ptsPermuted3D.col(posibleExit3Faces[iF][1])
					            +  barys(2) * ptsPermuted3D.col(posibleExit3Faces[iF][2]);

					/*Vec2 sum = (ptsProj2D.col(posibleExit3Faces[iF][0]) * barys[0]
						+ ptsProj2D.col(posibleExit3Faces[iF][1]) * barys[1]
						+ ptsProj2D.col(posibleExit3Faces[iF][2]) * barys[2]);*/

					// assert((ptsProj2D.col(posibleExit3Faces[iF][0]) * barys[0]
					// 	+ ptsProj2D.col(posibleExit3Faces[iF][1]) * barys[1]
					// 	+ ptsProj2D.col(posibleExit3Faces[iF][2]) * barys[2]).norm() < 1e-6);
					// assert(barys[0] > -0.01 && barys[0] < 1.01);
					// assert(barys[1] > -0.01 && barys[1] < 1.01);
					// assert(barys[2] > -0.01 && barys[2] < 1.01 );
					// assert((intersectPt - rayOrigin).cross(rayDir).norm() < 1e-5);
					

					if ( (intersectPt - rayOrigin).squaredNorm() > maxTraversalDis * maxTraversalDis)
					{
						statistics.stopReason = TraverseStopReason::passedMaximumDis;
						// std::cout << "Traversed distance: " << (intersectPt - rayOrigin).norm() 
						// 	<<" has passed through maximum traverse dis: " << maxTraversalDis
						// 	<<" at step " << steps << " | traversed set size:" <<  traversedTets.size() << "!\n";
						return false;
					}
				}
			}
		}
		++statistics.numTetsTraversed;

	}
	
	statistics.stopReason = TraverseStopReason::emptyStack;
	return false;
}

bool GAIA::TetMeshFEM::tetrahedralTraverseToLoopLess(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
	int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics)
{
	statistics.stopReason = TraverseStopReason::querySuccess;
	if (startTetId == targetTetId)
	{
		return true;
	}
	// get the ray-aligned 2D coordinate system
	Eigen::Matrix<FloatingType, 3, 2> axes;
	CuMatrix::buildOrthonormalBasis(rayDir.data(), axes.col(0).data(), axes.col(1).data());

	assert(axes.col(0).dot(axes.col(1)) < 1e-6);
	assert(axes.col(0).dot(rayDir) < 1e-6);
	assert(axes.col(1).dot(rayDir) < 1e-6);

	Eigen::Matrix<FloatingType, 2, 3> axesT = axes.transpose();

	// initialize the 2D coordinate
	Eigen::Matrix<FloatingType, 2, 4> ptsProj2D;

	// rearrange the ording of vertices thus the vertex across the incoming face is the last one
	// and the first 3 vertices are ordered as tet4Faces[startFaceId]
	projectTo2DCoordinates(ptsProj2D, startTetId, startFaceId, axesT, rayOrigin);

	Eigen::Vector4i possibleExitFace;
	exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);

	// figure out the outcoming Face
	std::vector<int32_t> candidateExitFaces;
	std::vector<int32_t> candidateExitFacesCorrespondingTets;

	candidateExitFaces.reserve(SIZE_CANDIDATE_FACE_STACK * 2);
	candidateExitFacesCorrespondingTets.reserve(SIZE_CANDIDATE_FACE_STACK * 2);


	for (size_t i = 0; i < 3; i++)
	{
		if (possibleExitFace(i))
		{
			// the first 3 vertices are ordered as tet4Faces[startFaceId]
			// thus tet4Faces[startFaceId][i] is the exit face id
			int32_t exitFaceId = tet4Faces[startFaceId][i];

			candidateExitFaces.push_back(exitFaceId);
			candidateExitFacesCorrespondingTets.push_back(startTetId);
		}
	}

	statistics.numTetsTraversed = 0;
	while (!candidateExitFaces.empty())
	{
		int32_t intersectedFaceId = candidateExitFaces.back();
		candidateExitFaces.pop_back();
		int32_t previousTetId = candidateExitFacesCorrespondingTets.back();
		candidateExitFacesCorrespondingTets.pop_back();

		int32_t intersectedFaceVIds[3] = {
			tetVIds()(tet4Faces[intersectedFaceId][0], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][1], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][2], previousTetId)
		};

		int32_t currentTetId = tetsNeighborTets()(intersectedFaceId, previousTetId);

		//std::cout << "Previous tet VIds: " << tetVIds.col(previousTetId).transpose() << std::endl;

		if (currentTetId == targetTetId)
		{
			statistics.stopReason = TraverseStopReason::querySuccess;
			return true;
		}

		if (currentTetId == -1)
			// reached the boundary
		{
			statistics.stopReason = TraverseStopReason::reachedBoundary;
			return false;
		}


		int32_t newVertexId = tetsXorSums()(currentTetId)
			^ intersectedFaceVIds[0]
			^ intersectedFaceVIds[1]
			^ intersectedFaceVIds[2];

		// the intersectedFaceId is not the incoming face id of intersectedFaceId
		// we need to query the id of the incoming face
		int32_t incomingFaceIdCurTet = 3;
		for (size_t iV = 0; iV < 3; iV++)
		{
			if (newVertexId == tetVIds()(iV, currentTetId)) {
				incomingFaceIdCurTet = iV;

				break;
			}
		}
		//std::cout << "Current tet VIds: " << tetVIds.col(currentTetId).transpose() << std::endl;

		// rearrange the ording of vertices thus the vertex across the incoming face (incomingFaceIdCurTet) is the last one
		// and the first 3 vertices are ordered as tet4Faces[startFaceId]
		projectTo2DCoordinates(ptsProj2D, currentTetId, incomingFaceIdCurTet, axesT, rayOrigin);
		int  numNumIntersectingFaces = exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);

		numNumIntersectingFaces = checkExitFaceForward(rayDir, currentTetId, incomingFaceIdCurTet, possibleExitFace);
		for (size_t iF = 0; iF < 3; iF++)
		{
			if (possibleExitFace(iF))
				// add intersected faces to candidate intersection face stack
			{
				// the exit face id is the actual vertex id (0~3), not how they are arranged in ptsProj2D
				// exit face id can be deduced from tet4Faces[incomingFaceIdCurTet]
				int32_t exitFaceId = tet4Faces[incomingFaceIdCurTet][iF];
				candidateExitFaces.push_back(exitFaceId);
				candidateExitFacesCorrespondingTets.push_back(currentTetId);

				int32_t nextTetId = tetsNeighborTets()(exitFaceId, currentTetId);
				if (nextTetId == targetTetId)
				{
					statistics.stopReason = TraverseStopReason::querySuccess;
					return true;
				}

				if ((!(statistics.numTetsTraversed % passThroughCheckSteps) || statistics.numTetsTraversed == 2)
					&& maxTraversalDis > 0.f)
					// check pass through
				{
					Vec3 barys, intersectPt;

					//std::cout << "ptsProj2D: \n" << ptsProj2D << std::endl;

					originBarycentric2DTriangle(ptsProj2D.col(posibleExit3Faces[iF][0]),
						ptsProj2D.col(posibleExit3Faces[iF][1]),
						ptsProj2D.col(posibleExit3Faces[iF][2]),
						barys);
					intersectPt = Vec3::Zero();

					Eigen::Matrix<FloatingType, 3, 4> ptsPermuted3D;
					copyRepermutedVerts(ptsPermuted3D, currentTetId, incomingFaceIdCurTet);

					for (size_t iV = 0; iV < 3; iV++)
					{
						intersectPt += barys(iV) * ptsPermuted3D.col(posibleExit3Faces[iF][iV]);
					}

					Vec2 sum = (ptsProj2D.col(posibleExit3Faces[iF][0]) * barys[0]
						+ ptsProj2D.col(posibleExit3Faces[iF][1]) * barys[1]
						+ ptsProj2D.col(posibleExit3Faces[iF][2]) * barys[2]);

					if ((intersectPt - rayOrigin).squaredNorm() > maxTraversalDis * maxTraversalDis)
					{
						statistics.stopReason = TraverseStopReason::passedMaximumDis;
						//std::cout << "Passed through maximum traverse dis!\n";
						return false;
					}
				}
			}
		}
		++statistics.numTetsTraversed;
	} // while

	statistics.stopReason = TraverseStopReason::emptyStack;
	return false;
}

bool GAIA::TetMeshFEM::tetrahedralTraverseToDynamic(const Vec3& rayOrigin, const Vec3& rayDir, const FloatingType maxTraversalDis, int32_t startTetId, int32_t startFaceId,
	int32_t targetTetId, FloatingType rayTriIntersectionEpsilon, TraverseStatistics& statistics)
{
	statistics.stopReason = TraverseStopReason::querySuccess;
	if (startTetId == targetTetId)
	{
		return true;
	}
	// get the ray-aligned 2D coordinate system
	Eigen::Matrix<FloatingType, 3, 2> axes;
	CuMatrix::buildOrthonormalBasis(rayDir.data(), axes.col(0).data(), axes.col(1).data());

	assert(axes.col(0).dot(axes.col(1)) < 1e-6);
	assert(axes.col(0).dot(rayDir) < 1e-6);
	assert(axes.col(1).dot(rayDir) < 1e-6);

	Eigen::Matrix<FloatingType, 2, 3> axesT = axes.transpose();

	// initialize the 2D coordinate

	Eigen::Matrix<FloatingType, 2, 4> ptsProj2D;

	int count = 0;
	// rearrange the ording of vertices thus the vertex across the incoming face is the last one
	// and the first 3 vertices are ordered as tet4Faces[startFaceId]
	projectTo2DCoordinates(ptsProj2D, startTetId, startFaceId, axesT, rayOrigin);

	Eigen::Vector4i possibleExitFace;
	exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);

	// figure out the outcoming Face
	std::vector<int32_t> candidateExitFaces;
	std::vector<int32_t> candidateExitFacesCorrespondingTets;

	candidateExitFaces.reserve(SIZE_CANDIDATE_FACE_STACK * 2);
	candidateExitFacesCorrespondingTets.reserve(SIZE_CANDIDATE_FACE_STACK *2);

#if defined TRAVERESED_TETS_SET
	std::unordered_set<int32_t> traversedTets;
	traversedTets.insert(startTetId);
#elif defined TRAVERESED_TETS_CIRCULAR_ARRAY
	CircularArray<int32_t, SIZE_TRAVERSED_CIRCULAR_ARRAY> traversedTets;
	traversedTets.push_back(startTetId);
#else
	std::vector<int32_t> traversedTets;
	traversedTets.reserve(SIZE_TRAVERSED_LIST_STACK * 2);
	traversedTets.push_back(startTetId);
#endif

#ifdef OUTPUT_TRAVERSED_TETS 
	statistics.statistics.traversedTetsOutput.clear();
	statistics.statistics.traversedTetsOutput.push_back(startTetId);
#endif


	for (size_t i = 0; i < 3; i++)
	{
		if (possibleExitFace(i))
		{
			// the first 3 vertices are ordered as tet4Faces[startFaceId]
			// thus tet4Faces[startFaceId][i] is the exit face id
			int32_t exitFaceId = tet4Faces[startFaceId][i];

			candidateExitFaces.push_back(exitFaceId);
			candidateExitFacesCorrespondingTets.push_back(startTetId);
		}
	};
	statistics.numTetsTraversed = 0;
	while (!candidateExitFaces.empty())
	{
		int32_t intersectedFaceId = candidateExitFaces.back();
		candidateExitFaces.pop_back();
		int32_t previousTetId = candidateExitFacesCorrespondingTets.back();
		candidateExitFacesCorrespondingTets.pop_back();

		int32_t intersectedFaceVIds[3] = {
			tetVIds()(tet4Faces[intersectedFaceId][0], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][1], previousTetId),
			tetVIds()(tet4Faces[intersectedFaceId][2], previousTetId)
		};

		int32_t currentTetId = tetsNeighborTets()(intersectedFaceId, previousTetId);

		//std::cout << "Previous tet VIds: " << tetVIds.col(previousTetId).transpose() << std::endl;

#ifdef OUTPUT_TRAVERSED_TETS 
		statistics.traversedTetsOutput.push_back(currentTetId);
#endif

		if (currentTetId == targetTetId)
		{
			statistics.stopReason = TraverseStopReason::querySuccess;
			return true;
		}

		if (currentTetId == -1)
			// reached the boundary
		{
			statistics.stopReason = TraverseStopReason::reachedBoundary;
			return false;
		}


#if defined TRAVERESED_TETS_SET

		if (traversedTets.find(currentTetId) != traversedTets.end())
			// formed loop, cut this branch here
		{
			continue;
		}
		else
		{
			traversedTets.insert(currentTetId);
		}
#elif defined TRAVERESED_TETS_CIRCULAR_ARRAY
		if (traversedTets.hasBackward(currentTetId))
		{
			continue;
		}
#ifndef ONLY_RECORD_BRANCHING_TETS
		traversedTets.push_back(currentTetId);
#endif // ONLY_RECORD_BRANCHING_TETS
#else
		bool loopDetected = false;
		int numTraversedTets = traversedTets.size();
		for (int iTet = 0; iTet < numTraversedTets; iTet++)
		{
			if (traversedTets[iTet] == currentTetId) {
				loopDetected = true;
				break;
			}
		}
		if (loopDetected)
		{
			continue;
		}
		traversedTets.push_back(currentTetId);
#endif

		int32_t newVertexId = tetsXorSums()(currentTetId)
			^ intersectedFaceVIds[0]
			^ intersectedFaceVIds[1]
			^ intersectedFaceVIds[2];

		// the intersectedFaceId is not the incoming face id of intersectedFaceId
		// we need to query the id of the incoming face
		int32_t incomingFaceIdCurTet = 3;
		for (size_t iV = 0; iV < 3; iV++)
		{
			if (newVertexId == tetVIds()(iV, currentTetId)) {
				incomingFaceIdCurTet = iV;

				break;
			}
		}
		//std::cout << "Current tet VIds: " << tetVIds.col(currentTetId).transpose() << std::endl;

		// rearrange the ording of vertices thus the vertex across the incoming face (incomingFaceIdCurTet) is the last one
		// and the first 3 vertices are ordered as tet4Faces[startFaceId]
		projectTo2DCoordinates(ptsProj2D, currentTetId, incomingFaceIdCurTet, axesT, rayOrigin);
		int  numNumIntersectingFaces = exitFaceSelection(ptsProj2D, possibleExitFace, rayTriIntersectionEpsilon);
#ifdef ONLY_RECORD_BRANCHING_TETS
		if (numNumIntersectingFaces > 1)
		{
			traversedTets.push_back(currentTetId);
		}
#endif // ONLY_RECORD_BRANCHING_TETS
		for (size_t iF = 0; iF < 3; iF++)
		{
			if (possibleExitFace(iF))
				// add intersected faces to candidate intersection face stack
			{
				// the exit face id is the actual vertex id (0~3), not how they are arranged in ptsProj2D
				// exit face id can be deduced from tet4Faces[incomingFaceIdCurTet]
				int32_t exitFaceId = tet4Faces[incomingFaceIdCurTet][iF];
				candidateExitFaces.push_back(exitFaceId);
				candidateExitFacesCorrespondingTets.push_back(currentTetId);
				
				int32_t nextTetId = tetsNeighborTets()(exitFaceId, currentTetId);
				if (nextTetId == targetTetId)
				{
					statistics.stopReason = TraverseStopReason::querySuccess;
					return true;
				}

				if ( ( !(statistics.numTetsTraversed % passThroughCheckSteps) || statistics.numTetsTraversed == 2)
					&& maxTraversalDis > 0.f)
					// check pass through
				{
					Vec3 barys, intersectPt;

					//std::cout << "ptsProj2D: \n" << ptsProj2D << std::endl;

					originBarycentric2DTriangle(ptsProj2D.col(posibleExit3Faces[iF][0]),
						ptsProj2D.col(posibleExit3Faces[iF][1]),
						ptsProj2D.col(posibleExit3Faces[iF][2]),
						barys);
					intersectPt = Vec3::Zero();

					Eigen::Matrix<FloatingType, 3, 4> ptsPermuted3D;
					copyRepermutedVerts(ptsPermuted3D, currentTetId, incomingFaceIdCurTet);

					for (size_t iV = 0; iV < 3; iV++)
					{
						intersectPt += barys(iV) * ptsPermuted3D.col(posibleExit3Faces[iF][iV]);
					}

					Vec2 sum = (ptsProj2D.col(posibleExit3Faces[iF][0]) * barys[0]
						+ ptsProj2D.col(posibleExit3Faces[iF][1]) * barys[1]
						+ ptsProj2D.col(posibleExit3Faces[iF][2]) * barys[2]);

					if ((intersectPt - rayOrigin).squaredNorm() > maxTraversalDis * maxTraversalDis)
					{
						statistics.stopReason = TraverseStopReason::passedMaximumDis;
						//std::cout << "Passed through maximum traverse dis!\n";
						return false;
					}
				}
			}
		}
		++statistics.numTetsTraversed;
	} // while

	statistics.stopReason = TraverseStopReason::emptyStack;
	return false;
}


void GAIA::TetMeshFEM::computeTopology(TetMeshMF* pTM_MF) 
{
	// scope of topologyLockGuard lock
	bool alreadyComputed = false;
	{
		std::lock_guard<std::mutex> topologyLockGuard(topologies_lock);
		std::string modelPath = pObjectParams->path;

		auto pTopoItem = topologies.find(modelPath);
		if (pTopoItem == topologies.end())
		{
			pTopology = std::make_shared<TetMeshTopology>();
			topologies.insert({ modelPath, pTopology });
			alreadyComputed = false;
		}
		else
		{
			pTopology = pTopoItem->second;
			alreadyComputed = true;
		}
	}

	if (!alreadyComputed)
	{
		pTopology->initialize(pTM_MF, pObjectParams);
	}
}
