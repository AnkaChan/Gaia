#pragma once
#include "VBDClothPhysics.h"

namespace GAIA {
	inline void setFixedPointCloth(VBDClothSimulationFramework& physics, const std::vector<int>& selectedMeshes, const std::vector<std::vector<int>>& selectedVerts)
	{
		for (size_t i = 0; i < selectedMeshes.size(); ++i)
		{
			auto iMesh = selectedMeshes[i];
			TriMeshFEM* pMesh = physics.baseTriMeshesForSimulation[iMesh].get();
			for (int vertexId : selectedVerts[i]) {
				Vec3 p = pMesh->vertex(vertexId);
				pMesh->fixedMask[vertexId] = true;

				if (std::find(pMesh->pObjectParams->fixedPoints.cbegin(), pMesh->pObjectParams->fixedPoints.cend(), vertexId)
					== pMesh->pObjectParams->fixedPoints.end())
				{
					pMesh->pObjectParams->fixedPoints.push_back(vertexId);
				}
			}

			/*if (physics.physicsParams().useGPU)
			{
				pTM->pTetMeshSharedBase->vertexFixedMaskBuffer->toGPU();
			}*/
		}
	}

	struct VBDClothBaseDeformer
	{
	public:
		VBDClothBaseDeformer() {};
		~VBDClothBaseDeformer() {};

		virtual void operator()(VBDClothSimulationFramework& physics, FloatingType curTime, int iFrame,
			int iSubstep, int iIter, FloatingType dt) = 0;

	private:

	};
	typedef std::shared_ptr<VBDClothBaseDeformer> ClothDeformerPtr;

	struct ClothDeformerRotator : public VBDClothBaseDeformer {
		ClothDeformerRotator(const Vec3 inO, const Vec3 inAxis, double inAngularVelocity, VBDClothSimulationFramework& physics,
			const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inRotationEndTime = -1) :
			o(inO), axis(inAxis), angularVelocity(inAngularVelocity), selectedMeshes(inSelectedMeshes),
			selectedVerts(inSelectedVerts), rotationEnd(inRotationEndTime)
		{
			setFixedPointCloth(physics, selectedMeshes, selectedVerts);

			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				TriMeshFEM* pMesh = physics.baseTriMeshesForSimulation[iMesh].get();
				roots.emplace_back();
				RootPs.emplace_back();
				for (int vertexId : selectedVerts[i]) {
					Vec3 p = pMesh->vertex(vertexId);

					Vec3 P = { p[0], p[1], p[2] };
					Vec3 OP = P - o;
					Vec3 root = OP.dot(axis) * axis;
					Vec3 RootP = P - root;

					roots.back().push_back(root);
					RootPs.back().push_back(RootP);
				}
			}
		}

		virtual void operator()(VBDClothSimulationFramework& physics, FloatingType curTime, int iFrame,
			int iSubstep, int iIter, FloatingType dt) {

			if (curTime >= rotationEnd && rotationEnd > 0)
			{
				curTime = rotationEnd;
			}
			double theta = curTime * angularVelocity;


			double ux = axis[0];
			double uy = axis[1];
			double uz = axis[2];

			Eigen::Matrix3f R;
			R << cos(theta) + ux * ux * (1 - cos(theta)), ux* uy* (1 - cos(theta)) - uz * sin(theta), ux* uz* (1 - cos(theta)) + uy * sin(theta),
				uy* ux* (1 - cos(theta)) + uz * sin(theta), cos(theta) + uy * uy * (1 - cos(theta)), uy* uz* (1 - cos(theta)) - ux * sin(theta),
				uz* ux* (1 - cos(theta)) - uy * sin(theta), uz* uy* (1 - cos(theta)) + ux * sin(theta), cos(theta) + uz * uz * (1 - cos(theta));

			// Eigen::Vector3d O = { o[0], o[1], o[2] };
			// Eigen::Vector3d N = { axis[0], axis[1], axis[2] };

			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				TriMeshFEM* pMesh = physics.baseTriMeshesForSimulation[iMesh].get();
				int iV = 0;
				for (int vId : selectedVerts[i]) {
					Vec3Block p = pMesh->vertex(vId);

					const Vec3& root = roots[i][iV];

					const Vec3& RootP = RootPs[i][iV];

					const Vec3 RootPRot = R * RootP;
					const Vec3 PRot = root + RootPRot;

					// Eigen::Vector3d PRot = root + RootP;
					//std::cout << "Rotate from: " << p.transpose() << " to " << PRot.transpose() << std::endl;

					p[0] = PRot[0];
					p[1] = PRot[1];
					p[2] = PRot[2];
					++iV;
				}
			}

		};

		Vec3 o;
		Vec3 axis;

		double angularVelocity;

		double rotationEnd;

		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		std::vector<std::vector<Vec3>> roots;
		std::vector<std::vector<Vec3>> RootPs;
	};

	struct ClothDeformerTranslater : public VBDClothBaseDeformer {
		ClothDeformerTranslater(const Vec3 inSpeed, VBDClothSimulationFramework& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, FloatingType inDeformationStartTime = 0.f, FloatingType inDeformationEndTime = -1) :
			speed(inSpeed), selectedMeshes(inSelectedMeshes), selectedVerts(inSelectedVerts), deformationStartTime(inDeformationStartTime), deformationEndTime(inDeformationEndTime)
		{
			setFixedPointCloth(physics, selectedMeshes, selectedVerts);
			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				TriMeshFEM* pMesh = physics.baseTriMeshesForSimulation[iMesh].get();
				originalPos.emplace_back();
				for (int vertexId : selectedVerts[i]) {
					Vec3 p = pMesh->vertex(vertexId);

					originalPos.back().push_back(p);
				}
			}
		}

		virtual void operator()(VBDClothSimulationFramework& physics, FloatingType curTime, int iFrame, int iSubstep, int iIter, FloatingType dt) {
			//curTime = curTime - deformationStartTime;
			//if (curTime > 0)
			//{

			//	if (curTime >= deformationEndTime - deformationStartTime && deformationEndTime > 0)
			//	{
			//		curTime = deformationEndTime - deformationStartTime;
			//	}


			//	for (size_t i = 0; i < selectedMeshes.size(); ++i)
			//	{
			//		auto iMesh = selectedMeshes[i];
			//		TriMeshFEM* pM = physics.baseTriMeshesForSimulation[iMesh].get();
			//		int iV = 0;
			//		for (int vId : selectedVerts[i]) {
			//			Vec3Block p = pM->vertex(vId);

			//			p = speed * curTime + originalPos[i][iV];;
			//			++iV;
			//		}
			//	}
			//}
			if (curTime >= deformationStartTime && (curTime < deformationEndTime || deformationEndTime < 0))
			{
				for (size_t i = 0; i < selectedMeshes.size(); ++i)
				{
					auto iMesh = selectedMeshes[i];
					TriMeshFEM* pM = physics.baseTriMeshesForSimulation[iMesh].get();
					int iV = 0;
					for (int vId : selectedVerts[i]) {
						Vec3Block p = pM->vertex(vId);
						p += speed * dt;
						++iV;
					}
				}
			}
		};

		Vec3 speed;
		FloatingType deformationStartTime;
		FloatingType deformationEndTime;
		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		std::vector<std::vector<Vec3>> originalPos;
	};


	inline ClothDeformerPtr loadClothDeformers(VBDClothSimulationFramework& physics, nlohmann::json& deformerParams)
	{
		std::string DeformerName;
		EXTRACT_FROM_JSON(deformerParams, DeformerName);

		if (DeformerName == "Rotator")
		{
			Vec3 center;
			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(deformerParams, center);
			Vec3 axis;
			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(deformerParams, axis);
			double angularVelocity;
			EXTRACT_FROM_JSON(deformerParams, angularVelocity);
			std::vector<int> selectedMeshes;
			EXTRACT_FROM_JSON(deformerParams, selectedMeshes);
			std::vector<std::vector<int>> selectedVertices;
			EXTRACT_FROM_JSON(deformerParams, selectedVertices);
			double rotationEndTime = -1;
			EXTRACT_FROM_JSON(deformerParams, rotationEndTime);

			return std::make_shared<ClothDeformerRotator>(center, axis, angularVelocity, physics, selectedMeshes, selectedVertices, rotationEndTime);
		}
		else if (DeformerName == "Translator") {
			Vec3 speed;
			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(deformerParams, speed);
			FloatingType deformationEndTime = -1;
			EXTRACT_FROM_JSON(deformerParams, deformationEndTime);
			FloatingType deformationStartTime = 0;
			EXTRACT_FROM_JSON(deformerParams, deformationStartTime);
			std::vector<int> selectedMeshes;
			EXTRACT_FROM_JSON(deformerParams, selectedMeshes);
			std::vector<std::vector<int>> selectedVertices;
			EXTRACT_FROM_JSON(deformerParams, selectedVertices);
			return std::make_shared<ClothDeformerTranslater>(speed, physics, selectedMeshes, selectedVertices, deformationStartTime, deformationEndTime);
		}
		else
		{
			std::cout << "Unrecognized deformer: " << DeformerName << "\n";
			assert(false);
			return nullptr;
		}

	}


}