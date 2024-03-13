#pragma once


#include "../Types/Types.h"
#include "VBD_BaseMaterial.h"
#include <memory>

#include "../ImplicitEuler/HessianAssebler.h"

#include "../TetMesh/TetMeshFEMShared.h"
#include "VBD_NeoHookeanGPU.h"



namespace GAIA {
	struct VBDTetMeshNeoHookean;

	struct ObjectParametersVBDNeoHookean : public ObjectParamsVBD
	{
		typedef std::shared_ptr<ObjectParametersVBDNeoHookean> SharedPtr;
		typedef ObjectParametersVBDNeoHookean* Ptr;

		ObjectParametersVBDNeoHookean() {
			materialType = NeoHookean;
		}
		FloatingType lmbd = 1e4;
		FloatingType miu = 1e4;
		FloatingType damping = 0;
		FloatingType dampingShear = 0;
		FloatingType dampingVolume = 0;

		virtual bool fromJson(nlohmann::json& objectJsonParams);
		virtual bool toJson(nlohmann::json& objectJsonParams);

	};

	struct VBDTetMeshNeoHookeanShared : public VBDBaseTetMeshShared
	{
		virtual void initializeGPUTetMesh(VBDTetMeshNeoHookean* pTetMesh, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, bool withAccelerator);
		virtual void initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);

	};

	class VBDTetMeshNeoHookean : public VBDBaseTetMesh {
	public:
		typedef std::shared_ptr<VBDTetMeshNeoHookean> SharedPtr;
		typedef VBDTetMeshNeoHookean* Ptr;

		virtual void solverIteration();
		virtual void evaluateInternalForce();

		virtual void VBDStep(int iV);
		FloatingType evaluateVertexMeritEnergy(int iV, FloatingType& meInertia, FloatingType& meElastic_bending);
		virtual void accumlateMaterialForceAndHessian(int iV, Vec3& force, Mat3& hessian);
		virtual void accumlateMaterialForceAndHessian2(int iV, Vec3& force, Mat3& hessian);
		virtual void accumlateMaterialForce(int iV, Vec3& force);
		virtual void computeElasticForceHessian(int tetId, FloatingType& energy, Vec12& force, Mat12& hessian);
		virtual void computeElasticForceHessianDouble(int tetId, double& energy, Eigen::Vector<double, 12>& force, Eigen::Matrix<double, 12, 12>& hessian);
		virtual void computeElasticGradientHessianF(const Mat3& F, FloatingType& energy, Vec9& gradient, Mat9& hessian);
		virtual void computeElasticGradientHessianFDouble(const Eigen::Matrix3d& F, double& energy, Eigen::Vector<double, 9>& gradient, Eigen::Matrix<double, 9, 9>& hessian);
		void validateElasticForceHessian(int tetId);
		void validateElasticGradientHessianF(int tetId);

		template<typename T>
		inline void computeElasticEnergy(int tetId, T& energy)
		{
			const auto& material = ObjectParametersMaterial();
			T miu = material.miu;
			T lmbd = material.lmbd;
			T a = 1 + miu / lmbd;
			T A = tetRestVolume(tetId);

			Eigen::Matrix3<T> DmInv = getDmInv(tetId).cast<T>();
			Mat3 DsFloat;
			computeDs(DsFloat, tetId);
			Eigen::Matrix3<T> Ds = DsFloat.cast<T>();
			Eigen::Matrix3<T> F = Ds * DmInv;
			const T detF = F.determinant();
			const T k = detF - a;
			energy = 0.5 * A * (lmbd * k * k + miu * (F.squaredNorm() - 3));

			T restEnergy = 0.5 * A * lmbd * SQR(1.f - a);

			energy -= restEnergy;
		}

		const ObjectParametersVBDNeoHookean& ObjectParametersMaterial();

		FloatingType backTracingLineSearchVBD(IdType vId, const Vec3& dx, FloatingType E0, FloatingType alpha,
			FloatingType c, FloatingType tau, int maxNumIters);

		virtual void initializeGPUMesh();
		virtual VBDBaseTetMeshGPU* getGPUMesh();
		virtual VBDBaseTetMeshGPU* getGPUMesh_forCPUDebug();

		virtual void syncToGPU(bool sync = false, cudaStream_t stream = 0);
		virtual void syncToCPU(bool sync = false, cudaStream_t stream = 0);
		virtual void syncToGPUVertPosOnly(bool sync = false, cudaStream_t stream = 0);
		virtual void syncToCPUVertPosOnly(bool sync = false, cudaStream_t stream = 0);
		virtual void setGPUMeshActiveness(bool activeForCollision, bool activeForMaterialSolve);


	public:
		std::shared_ptr<VBDTetMeshNeoHookeanShared> pTetMeshShared;
		std::shared_ptr<DeviceClassBuffer<VBDTetMeshNeoHookeanGPU>> pTetMeshGPUBuffer;
		std::shared_ptr<VBDTetMeshNeoHookeanGPU> pTetMeshGPUBufferOnCPU;
		std::shared_ptr<VBDTetMeshNeoHookeanGPU> pTetMeshGPUBuffer_forCPUDebug;

	};

	inline const ObjectParametersVBDNeoHookean& GAIA::VBDTetMeshNeoHookean::ObjectParametersMaterial()
	{
		return  *(ObjectParametersVBDNeoHookean*)(pObjParamsVBD.get());
	}

	inline FloatingType VBDTetMeshNeoHookean::backTracingLineSearchVBD(IdType vId, const Vec3& dx, FloatingType E0,
		FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters)
	{
		FloatingType m = dx.squaredNorm();

		const Vec3 orgPos = vertex(vId);

		FloatingType bestAlpha = 0.f;
		FloatingType bestEnergy = E0;

		FloatingType meInertia = 0;
		FloatingType meElastic_elasitic = 0;

		for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
		{
			vertex(vId) = alpha * dx + orgPos;

			FloatingType e = evaluateVertexMeritEnergy(vId, meInertia, meElastic_elasitic);

			if (e < bestEnergy)
			{
				bestAlpha = alpha;
				bestEnergy = e;
			}

			// first Wolfe condition 
			if (e < E0 - alpha * c * m)
			{
				// std::cout << "step size for vertex " << vId << ": " << alpha << "\n";
				break;
			}
			else
			{
				alpha = alpha * tau;
			}
		}

		return bestEnergy;
	}





}