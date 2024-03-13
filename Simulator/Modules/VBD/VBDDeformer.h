#pragma once 
#include "VBDPhysics.h"

namespace GAIA {

	struct VBDBaseDeformer
	{
	public:
		VBDBaseDeformer() {};
		~VBDBaseDeformer() {};

		virtual void operator()(VBDPhysics& physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) = 0;

	private:

	};

	typedef std::shared_ptr<VBDBaseDeformer> DeformerPtr;

	inline void setFixed(VBDPhysics& physics, const std::vector<int>& selectedMeshes, const std::vector<std::vector<int>>& selectedVerts) 
	{
		for (size_t i = 0; i < selectedMeshes.size(); ++i)
		{
			auto iMesh = selectedMeshes[i];
			VBDBaseTetMesh* pTM = physics.tMeshes[iMesh].get();
			for (int vertexId : selectedVerts[i]) {
				Vec3 p = pTM->vertex(vertexId);
				pTM->verticesCollisionDetectionEnabled[vertexId] = false;
				pTM->fixedMask[vertexId] = true;

				if (std::find(pTM->pObjectParams->fixedPoints.cbegin(), pTM->pObjectParams->fixedPoints.cend(), vertexId)
					== pTM->pObjectParams->fixedPoints.end())
				{
					pTM->pObjectParams->fixedPoints.push_back(vertexId);
				}
			}

			if (physics.physicsParams().useGPU)
			{
				pTM->pTetMeshSharedBase->vertexFixedMaskBuffer->toGPU();
			}
		}
	}

	struct DeformerRotator : public VBDBaseDeformer {
		DeformerRotator(const Vec3 inO, const Vec3 inAxis, double inAngularVelocity, VBDPhysics& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inRotationEndTime = -1) :
			o(inO), axis(inAxis), angularVelocity(inAngularVelocity), selectedMeshes(inSelectedMeshes),
			selectedVerts(inSelectedVerts), rotationEnd(inRotationEndTime)
		{
			setFixed(physics, selectedMeshes, selectedVerts);

			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				VBDBaseTetMesh* pTM = physics.tMeshes[iMesh].get();
				roots.emplace_back();
				RootPs.emplace_back();
				for (int vertexId : selectedVerts[i]) {
					Vec3 p = pTM->vertex(vertexId);

					Vec3 P = { p[0], p[1], p[2] };
					Vec3 OP = P - o;
					Vec3 root = OP.dot(axis) * axis;
					Vec3 RootP = P - root;

					roots.back().push_back(root);
					RootPs.back().push_back(RootP);
				}
			}
		}

		virtual void operator()(VBDPhysics& physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {

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
				TetMeshFEM* pTM = physics.tMeshes[iMesh].get();
				int iV = 0;
				for (int vId : selectedVerts[i]) {
					Vec3Block p = pTM->vertex(vId);

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


	struct DeformerTranslater : public VBDBaseDeformer {
		DeformerTranslater(const Vec3 inSpeed, VBDPhysics& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inDeformationEndTime = -1) :
			speed(inSpeed), selectedMeshes(inSelectedMeshes), selectedVerts(inSelectedVerts), deformationEndTime(inDeformationEndTime)
		{
			setFixed(physics, selectedMeshes, selectedVerts);
			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				VBDBaseTetMesh* pTM = physics.tMeshes[iMesh].get();
				originalPos.emplace_back();
				for (int vertexId : selectedVerts[i]) {
					Vec3 p = pTM->vertex(vertexId);

					originalPos.back().push_back(p);
				}
			}
		}

		virtual void operator()(VBDPhysics& physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
			if (curTime >= deformationEndTime && deformationEndTime > 0)
			{
				curTime = deformationEndTime;
			}
			for (size_t i = 0; i < selectedMeshes.size(); ++i)
			{
				auto iMesh = selectedMeshes[i];
				TetMeshFEM* pTM = physics.tMeshes[iMesh].get();
				int iV = 0;
				for (int vId : selectedVerts[i]) {
					Vec3Block p = pTM->vertex(vId);

					p = speed * curTime + originalPos[i][iV];;
					++iV;
				}
			}

		};

		Vec3 speed;
		FloatingType deformationEndTime ;
		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		std::vector<std::vector<Vec3>> originalPos;
	};


	inline DeformerPtr loadDeformers(VBDPhysics& physics, nlohmann::json& deformerParams)
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

			return std::make_shared<DeformerRotator>(center, axis, angularVelocity, physics, selectedMeshes, selectedVertices, rotationEndTime);
		}
		else if (DeformerName == "Translator") {
			Vec3 speed;
			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(deformerParams, speed);
			double deformationEndTime = -1;
			EXTRACT_FROM_JSON(deformerParams, deformationEndTime);
			std::vector<int> selectedMeshes;
			EXTRACT_FROM_JSON(deformerParams, selectedMeshes);
			std::vector<std::vector<int>> selectedVertices;
			EXTRACT_FROM_JSON(deformerParams, selectedVertices);
			return std::make_shared<DeformerTranslater>(speed, physics, selectedMeshes, selectedVertices, deformationEndTime);
		}
		else
		{
			std::cout << "Unrecognized deformer: " << DeformerName << "\n";
			assert(false);
			return nullptr;
		}

	}
}