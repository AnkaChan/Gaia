#include "FileIO.h"
#include "../Framework/BasePhysicsFramework.h"
#include <MeshFrame/Utility/IO.h>
#include <MeshFrame/Utility/Str.h>

using namespace GAIA;

void GAIA::saveAsPly(const std::string filePath, const std::vector<std::array<FloatingType, 3>>& verts, const std::vector<std::array<int, 3>>& faces) {
	//std::ofstream ofsPLY(filePath);
	//std::stringstream ss;

	FILE* fp = fopen(filePath.c_str(), "w");

	fprintf(fp, "ply\nformat ascii 1.0\ncomment Mocap generated\nelement vertex %d\n", verts.size());
	fprintf(fp, "property float x\nproperty float y\nproperty float z\n");

	if (faces.size() != 0) {

		fprintf(fp, "element face %d \nproperty list uchar int vertex_indices\n", faces.size());
	}
	fprintf(fp, "end_header\n");

	//ss << ;
	//fprintf(fp, "%s", ss.str());

	for (int i = 0; i < verts.size(); i++)
	{
		// ss << verts[i][0] << " " << verts[i][1] << " " << verts[i][2] << "\n";
		fprintf(fp, "%f %f %f\n", verts[i][0], verts[i][1], verts[i][2]);
	}

	for (int i = 0; i < faces.size(); i++) {
		//ss << faces[i].size();
		//for (size_t j = 0; j < faces[i].size(); j++)
		//{
		//	ss << " " << faces[i][j];
		//}
		//ss << "\n";
		fprintf(fp, "3 %d %d %d\n", faces[i][0], faces[i][1], faces[i][2]);

	}


	////fputs("This is testing for fputs...\n", fp);
	fclose(fp);

	//ofsPLY << ss.str();

	//ofsPLY.close();
}

void GAIA::writeAllToBinary(const char* output, std::vector<TetMeshFEM::SharedPtr>& tetMeshes)
{
	if (!tetMeshes.size())
	{
		return;
	}
	int mAllVerts = 0;
	std::vector<FloatingType> verts;

	verts.reserve(3 * tetMeshes.size() * tetMeshes[0]->numVertices());

	for (int iMesh = 0; iMesh < tetMeshes.size(); ++iMesh) {
		TetMeshFEM::SharedPtr& pTMesh = tetMeshes[iMesh];

		for (int iV = 0; iV < pTMesh->surfaceVIds().size(); ++iV)
		{
			auto v = pTMesh->surfaceVertex(iV);

			verts.push_back(v[0]);
			verts.push_back(v[1]);
			verts.push_back(v[2]);
		}
	}

	std::ofstream out;
	out.open(output, std::ios::out | std::ios::binary);
	for (int iV = 0; iV < verts.size(); iV++)
	{
		out.write(reinterpret_cast<char*>(&verts[iV]), sizeof(float));
	}
	out.close();
}

void GAIA::writeAllToPLY(const char* output, std::vector<TetMeshFEM::SharedPtr>& tetMeshes, bool saveAllInOneFile, bool addModelName)
{
	if (!saveAllInOneFile)
	{
		for (int iMesh = 0; iMesh < tetMeshes.size(); ++iMesh) {
			TetMeshFEM::SharedPtr& pTMesh = tetMeshes[iMesh];
			MF::IO::FileParts fp(output);

			std::ostringstream aSs;
			aSs << std::setfill('0') << std::setw(6) << iMesh;
			std::string modelNumber = aSs.str();

			std::string outFileCurModel = fp.path + "/" + fp.name + "_Model_" + modelNumber;
			if (addModelName)
			{
				std::string modelPath = pTMesh->pObjectParams->path;
				MF::IO::FileParts fp_modelPath(modelPath);
				outFileCurModel = outFileCurModel + "_" + fp_modelPath.name;
			}
			outFileCurModel = outFileCurModel + fp.ext;

			pTMesh->saveAsPLY(outFileCurModel.c_str());
		}
	}
	else {
		int mAllVerts = 0;
		std::vector<std::array<FloatingType, 3>> verts;
		std::vector<std::array<int, 3>> faces;

		for (int iMesh = 0; iMesh < tetMeshes.size(); ++iMesh) {
			TetMeshFEM::SharedPtr& pTMesh = tetMeshes[iMesh];

			for (int iV = 0; iV < pTMesh->surfaceVIds().size(); ++iV)
			{
				auto v = pTMesh->surfaceVertex(iV);
				std::array<FloatingType, 3> pt = { v[0], v[1], v[2] };

				verts.push_back(pt);
			}


			for (int iF = 0; iF < pTMesh->surfaceFacesSurfaceMeshVIds().cols(); ++iF)
			{
				std::array<int, 3> fVIds;
				for (size_t iFV = 0; iFV < 3; iFV++)
				{
					fVIds[iFV] = pTMesh->surfaceFacesSurfaceMeshVIds()(iFV, iF) + mAllVerts;
				}
				faces.push_back(fVIds);

			}

			mAllVerts += pTMesh->numSurfaceVerts();

		}
		//auto t1 = high_resolution_clock::now();

		// PathFinder::M mesh;
		// mesh.readVFList(&verts, &faces);
		// mesh.write_ply(output, 3);

		saveAsPly(output, verts, faces);
		//auto t2 = high_resolution_clock::now();
		//double timeConsumptionWritingPLY = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
		//std::cout << "timeConsumptionWritingPLY: " << timeConsumptionWritingPLY << std::endl;
	}
}

void GAIA::writeAllToPLY(const char* output, std::vector<TriMeshFEM::SharedPtr>& triMeshes, bool saveAllInOneFile)
{
	if (!saveAllInOneFile)
	{
		for (int iMesh = 0; iMesh < triMeshes.size(); ++iMesh) {
			TriMeshFEM::SharedPtr& pMesh = triMeshes[iMesh];
			MF::IO::FileParts fp(output);

			std::ostringstream aSs;
			aSs << std::setfill('0') << std::setw(4) << iMesh;
			std::string modelNumber = aSs.str();

			std::string outFileCurModel = fp.path + "/" + "Model_" + modelNumber ;
		
			outFileCurModel = outFileCurModel + "_" + fp.name + fp.ext;

			pMesh->saveAsPLY(outFileCurModel.c_str());
		}
	}
	else {
		int mAllVerts = 0;
		std::vector<std::array<FloatingType, 3>> verts;
		std::vector<std::array<int, 3>> faces;

		for (int iMesh = 0; iMesh < triMeshes.size(); ++iMesh) {
			TriMeshFEM::SharedPtr& pMesh = triMeshes[iMesh];

			for (int iV = 0; iV < pMesh->numVertices(); ++iV)
			{
				auto v = pMesh->positions().col(iV);
				std::array<FloatingType, 3> pt = { v[0], v[1], v[2] };

				verts.push_back(pt);
			}


			for (int iF = 0; iF < pMesh->numFaces(); ++iF)
			{
				std::array<int, 3> fVIds;
				for (size_t iFV = 0; iFV < 3; iFV++)
				{
					fVIds[iFV] = pMesh->facePosVId(iF, iFV) + mAllVerts;
				}
				faces.push_back(fVIds);

			}

			mAllVerts += pMesh->numVertices();

		}
		//auto t1 = high_resolution_clock::now();

		// PathFinder::M mesh;
		// mesh.readVFList(&verts, &faces);
		// mesh.write_ply(output, 3);

		saveAsPly(output, verts, faces);
		//auto t2 = high_resolution_clock::now();
		//double timeConsumptionWritingPLY = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
		//std::cout << "timeConsumptionWritingPLY: " << timeConsumptionWritingPLY << std::endl;
	}
}
