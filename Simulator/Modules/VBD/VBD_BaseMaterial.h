#pragma once
#include "../TetMesh/TetMeshFEM.h"
#include "VBDPhysicsParameters.h"

#include "../Modules/Utility/Logger.h"
#include "VBD_BaseMeshGPU.h"
#include <math.h>

#include "../TetMesh/TetMeshFEMShared.h"

#define BLOCK_DIM3x1_AT(v, iV) (v.block<POINT_VEC_DIMS, 1>(POINT_VEC_DIMS * iV, 0))

//#define OUTPUT_INITIALIZATION_GRAV_NORM

namespace GAIA {

	struct VBDPhysics;
	struct BasePhysicFramework;
	class VBDBaseTetMesh;

	struct ObjectParamsVBD : ObjectParams
	{
		typedef std::shared_ptr<ObjectParamsVBD> BaseSharedPtr;
		typedef ObjectParamsVBD* BasePtr;

		typedef std::shared_ptr<ObjectParamsVBD> SharedPtr;
		typedef ObjectParamsVBD* Ptr;

		FloatingType frictionDynamic = 0.1f;
		FloatingType frictionEpsV = 1.0f;

		FloatingType exponentialVelDamping = 1.0f;
		FloatingType constantVelDamping = 0.0f;

		bool localDescent = true;

		// 0: inertia, 1: previous position; 
		// 2: inertia without external force, 
		// 3: half inertia
		// 4: initRatio_g of inertia
		// 5: symplectic Euler
		// 6: adaptive
		int initializationType = 0;

		FloatingType initRatio_g = 0.5f;

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};


	struct VBDBaseTetMeshShared : public TetMeshFEMShared {
		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetForceAndHessiansBuffer;
		ManagedBuffer<int8_t>::SharedPtr activeCollisionMaskBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr inertiaBuffer;
#ifdef GPU_JACOBI_DX
		ManagedBuffer<FloatingTypeGPU>::SharedPtr positionsNewBuffer;
#endif // GPU_JACOBI_DX

		// for acceleration
		ManagedBuffer<FloatingTypeGPU>::SharedPtr positionsPrevIterBuffer;


		ManagedBuffer<CollisionDataGPU>::SharedPtr collisionDataBuffer;

		virtual void initializeGPUTetMesh(VBDBaseTetMesh* pTetMesh, VBDBaseTetMeshGPU* pTetMeshGPU, bool withAccelerator);
		virtual void initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);

		virtual void syncToGPU(bool sync = true, cudaStream_t stream = 0);
		virtual void syncToCPU(bool sync = true, cudaStream_t stream = 0);

		CollisionDataGPU& getCollisionDataCPUBuffer(int32_t iVert) { return collisionDataBuffer->getCPUBuffer()[iVert]; }
		CollisionRelationGPU& getCollisionRelationCPUBuffer(int32_t iVert, int32_t iRelation)
		{
			return getCollisionDataCPUBuffer(iVert).collisionRelations[iRelation];
		}

		int32_t& numCollisionRelation(int32_t iVert) { return collisionDataBuffer->getCPUBuffer()[iVert].numCollisionRelations; }

	};

	/*
	* A base material class for all other materials
	*
	*/
	class VBDBaseTetMesh : public TetMeshFEM {
	public:
		typedef std::shared_ptr<VBDBaseTetMesh> BaseSharedPtr;
		typedef VBDBaseTetMesh* BasePtr;

		typedef std::shared_ptr<VBDBaseTetMesh> SharedPtr;
		typedef VBDBaseTetMesh* Ptr;

		VBDBaseTetMesh()
		{
		};
		~VBDBaseTetMesh() {};

		virtual void computeInertia();
		virtual FloatingType computeInertiaEnergy();
		virtual void forwardStep();

		virtual void applyInitialStep();

		virtual void forwardStepSimplecticEuler();

		Eigen::Block<TVerticesMat, 3, -1> velocitiesPrev();

		virtual void solve();

		virtual void solverIteration() = 0;

		virtual void updateVelocity();

		virtual void handleForceConstraints();
		virtual void handleVelocityConstraint();
		virtual void evaluateExternalForce();;

		virtual void initializeGPUMesh() = 0;
		virtual VBDBaseTetMeshGPU* getGPUMesh() = 0;
		virtual VBDBaseTetMeshGPU* getGPUMesh_forCPUDebug() = 0;
		virtual void syncToGPU(bool sync = false, cudaStream_t stream = 0) = 0;
		virtual void syncToCPU(bool sync = false, cudaStream_t stream = 0) = 0;

		virtual void syncToGPUVertPosOnly(bool sync = false, cudaStream_t stream = 0) = 0;
		virtual void syncToCPUVertPosOnly(bool sync = false, cudaStream_t stream = 0) = 0;
		virtual void setGPUMeshActiveness(bool activeForCollision, bool activeForMaterialSolve) = 0;

		virtual void initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF,
			VBDPhysicsParameters::SharedPtr inPhysicsParams, BasePhysicFramework* in_pPhysicsFramework);

		virtual void evaluateInternalForce() = 0;
		// constant force is inrelavant to position and velocity
		virtual void evaluateConstantForce(Eigen::Ref<VecDynamic> internalForceOut);

		void accumlateInertiaForceAndHessian(int iV, Vec3& force, Mat3& hessian);

		template <typename Derived3x1, typename Derived3x3_1, typename Derived3x3_2>
		void computeDeformationGradient(const Eigen::MatrixBase<Derived3x1>& p1, const Eigen::MatrixBase<Derived3x1>& p2, const Eigen::MatrixBase<Derived3x1>& p3,
			const Eigen::MatrixBase<Derived3x1>& p4, const Eigen::MatrixBase<Derived3x3_1>& DsInv, Eigen::MatrixBase<Derived3x3_2>& F);

		FloatingType evaluatedMeritEnergy(FloatingType& meInertia, FloatingType& meElastic_elastic);
		FloatingType evaluatedAvgGradNorm();

		/*
		* Using the chain rule, given a function f's gradient of flatten F, we convert it to the gradient of flatten X
		* F is flattened in a row major manner
		* i.e., F_flatten = transpose([F(1,1) F(1,2) F(1,3) F(2,1) F(2,2) F(2,3) F(3,1) F(3,2) F(3,3)]);
		*/
		void chainRuleFGrad_ToXGrad(const Mat3& DsInv, const Vec9& FGrad, Vec3& p1Grad, Vec3& p2Grad, Vec3& p3Grad, Vec3& p4Grad);

		virtual Eigen::Block<VecDynamic> getPositionsBlock(VecDynamic& concatenatedVec) {
			return concatenatedVec.block(0, 0, 3 * numVertices(), 1);
		}
		virtual Eigen::Block<VecDynamic> getVelocityBlock(VecDynamic& concatenatedVec) {
			return concatenatedVec.block(3 * numVertices(), 0, 3 * numVertices(), 1);
		};

		//protected:
	public:
		ObjectParamsVBD::SharedPtr pObjParamsVBD;

		// infos related to face-based energy
		// Mat12x tetInternalForces;
		// size: 3 * numFace(),
		// 3*iFace ~ 3*IFace + 2 entries store the diagonal blocks of the Hessian matrix of the energy wrt. the triangle's vertices

		TVerticesMat vertexInternalForces;
		TVerticesMat vertexExternalForces;
		TVerticesMat mVelocitiesPrev;
		TVerticesMat inertia;
		TVerticesMat gradient;

		TVerticesMat acceletration;

		// records which parallel Group each vertex belongs 
		VecDynamicI vertexParallelGroups;
		// size: numVertices(), 
		// iV-th element stores the diagonal blocks of the Hessian matrix of the energy wrt. the iV-th vertex

		VecDynamic meritEnergy_Inertia_PerVertex;
		VecDynamic meritEnergyElastic_PerTet;

		VecDynamicBool activeCollisionMask;
		// record whether the point has penetrated in the beginning of a frame
		// used to separate dcd/ccd scheme
		VecDynamicBool penetratedMask;

		VBDPhysicsParameters::SharedPtr pPhysicsParams;
		BasePhysicFramework* pPhysicsFramework;

		std::shared_ptr<VBDBaseTetMeshShared> pTetMeshSharedBase;

		bool hasVelocitiesPrev = false;
		bool hasApproxAcceleration = false;
	};

	typedef std::shared_ptr<VBDBaseTetMesh> BaseMaterialPtr;

	inline void VBDBaseTetMesh::computeInertia()
	{
		CFloatingType dt = pPhysicsParams->dt;
		evaluateExternalForce();
		//cpu_parallel_for(0, numVertices(), [&](int iV)
		for (int iV = 0; iV < numVertices(); iV++)
		{
			mVelocity.col(iV) += dt * vertexExternalForces.col(iV) * vertexInvMass(iV);
		}
		//);
		//velocities() += dt * vertexExternalForces * vertexInvMass;


		if (pPhysicsParams->associateGravityWithInertia)
		{
			mVelocity.colwise() += dt * pPhysicsParams->gravity;
		}

		handleVelocityConstraint();
		inertia = positions() + dt * velocities();
	}

	inline FloatingType VBDBaseTetMesh::computeInertiaEnergy()
	{
		//std::cout << inertia.rows() << " x " << inertia.cols() << " | "
		//	<< vertices().rows() << " x " << vertices().cols() << std::endl;
		//auto diff = inertia - vertices();
		//auto a = diff.colwise().squaredNorm().array();

		//std::cout << "a: " << a.rows() << " x " << a.cols() << " | "
		//	<< "vertexMass: " << vertexMass.array().rows() << " x " << vertexMass.array().cols() << std::endl;
		FloatingType energy = ((inertia - vertices()).colwise().squaredNorm().array() * vertexMass.transpose().array()).sum();
		//FloatingType energy = 0;
		//for (int iV = 0; iV < numVertices(); iV++)
		//{
		//	energy += vertexMass(iV) * (inertia.col(iV) - vertex(iV)).squaredNorm();
		//}
		energy *= 0.5 * pPhysicsParams->dtSqrReciprocal;
		return energy;
	}



	inline void VBDBaseTetMesh::applyInitialStep()
	{
		positionsPrev() = positions();
		if (hasVelocitiesPrev)
		{
			acceletration = (mVelocity - mVelocitiesPrev) / pPhysicsParams->dt;
			hasApproxAcceleration = true;
		}
		else
		{
			hasVelocitiesPrev = true;
		}
		velocitiesPrev() = velocities();
		computeInertia();
		forwardStep();
	}


	inline Eigen::Block<TVerticesMat, 3, -1> VBDBaseTetMesh::velocitiesPrev()
	{
		return mVelocitiesPrev.block<3, -1>(0, 0, 3, numVertices());;
	}

	inline void VBDBaseTetMesh::solve()
	{
		applyInitialStep();

		evaluateExternalForce();

		for (size_t iteration = 0; iteration < pPhysicsParams->iterations; iteration++)
		{
			solverIteration();

		}

		updateVelocity();
	}

	inline void VBDBaseTetMesh::updateVelocity()
	{
		velocities() = (positions() - positionsPrev()) / pPhysicsParams->dt;
		for (size_t iFixedPoint = 0; iFixedPoint < pObjParamsVBD->fixedPoints.size(); iFixedPoint++)
		{
			mVelocity.col(pObjParamsVBD->fixedPoints[iFixedPoint]).setZero();
		}
	}

	inline void VBDBaseTetMesh::handleForceConstraints()
	{
		for (size_t iFixedPoint = 0; iFixedPoint < pObjParamsVBD->fixedPoints.size(); iFixedPoint++)
		{
			vertexInternalForces.col(pObjParamsVBD->fixedPoints[iFixedPoint]).setZero();
			vertexExternalForces.col(pObjParamsVBD->fixedPoints[iFixedPoint]).setZero();
		}
	}

	inline void VBDBaseTetMesh::handleVelocityConstraint()
	{
		for (size_t iFixedPoint = 0; iFixedPoint < pObjParamsVBD->fixedPoints.size(); iFixedPoint++)
		{
			velocities().col(pObjParamsVBD->fixedPoints[iFixedPoint]).setZero();
		}
	}

	inline void VBDBaseTetMesh::evaluateExternalForce()
	{
		vertexExternalForces.setZero();
	}


	inline void VBDBaseTetMesh::initialize(ObjectParams::SharedPtr inMaterialParams,
		std::shared_ptr<TetMeshMF> pTM_MF, VBDPhysicsParameters::SharedPtr inPhysicsParams, BasePhysicFramework* in_pPhysicsFramework)
	{
		TetMeshFEM::initialize(inMaterialParams, pTM_MF);
		pPhysicsParams = inPhysicsParams;
		pPhysicsFramework = in_pPhysicsFramework;
		pObjParamsVBD = std::static_pointer_cast<ObjectParamsVBD>(inMaterialParams);

		vertexExternalForces.resize(3, numVertices());
		vertexInternalForces.resize(3, numVertices());

		inertia.resize(3, numVertices());

		meritEnergy_Inertia_PerVertex.resize(numVertices());
		meritEnergyElastic_PerTet.resize(numTets());

		gradient.resize(3, numVertices());
		gradient.setZero();

		mVelocitiesPrev.resize(3, numVertices() + 1);
		mVelocitiesPrev = mVelocity;

		activeCollisionMask.resize(numVertices());
		activeCollisionMask.setZero();

		penetratedMask.resize(numVertices());
		penetratedMask.setZero();

		vertexParallelGroups.resize(numVertices());

		if (pObjParamsVBD->initializationType == 6)
		{
			acceletration.resize(3, numVertices() + 1);
		}
	}

	inline void VBDBaseTetMesh::evaluateConstantForce(Eigen::Ref<VecDynamic> forceOut)
	{
		// handle gravity
		for (size_t iVert = 0; iVert < numVertices(); iVert++)
		{
			BLOCK_DIM3x1_AT(forceOut, iVert) = pPhysicsParams->gravity * vertexMass(iVert);
		}
	}

	inline void VBDBaseTetMesh::accumlateInertiaForceAndHessian(int iV, Vec3& force, Mat3& hessian)
	{
		CFloatingType dtSqrReciprocal = pPhysicsParams->dtSqrReciprocal;

		force += vertexMass(iV) * (inertia.col(iV) - vertex(iV)) * (dtSqrReciprocal);

		hessian += Mat3::Identity() * vertexMass(iV) * (dtSqrReciprocal);
	}

	inline FloatingType VBDBaseTetMesh::evaluatedMeritEnergy(FloatingType& meInertia, FloatingType& meElastic_elastic)
	{

		return FloatingType();
	}

	inline FloatingType VBDBaseTetMesh::evaluatedAvgGradNorm()
	{
		return FloatingType();
	}

	template <typename Derived3x1, typename Derived3x3_1, typename Derived3x3_2>
	void VBDBaseTetMesh::computeDeformationGradient(const Eigen::MatrixBase<Derived3x1>& p1, const Eigen::MatrixBase<Derived3x1>& p2, const Eigen::MatrixBase<Derived3x1>& p3,
		const Eigen::MatrixBase<Derived3x1>& p4, const Eigen::MatrixBase<Derived3x3_1>& DsInv, Eigen::MatrixBase<Derived3x3_2>& F)
	{
		Mat3 Dm;
		Dm.block<3, 1>(0, 0) = p2 - p1;
		Dm.block<3, 1>(0, 1) = p3 - p1;
		Dm.block<3, 1>(0, 2) = p4 - p1;

		F = Dm * DsInv;
	}


	inline void GAIA::VBDBaseTetMesh::chainRuleFGrad_ToXGrad(const Mat3& DsInv, const Vec9& FGrad, Vec3& p1Grad, Vec3& p2Grad, Vec3& p3Grad, Vec3& p4Grad)
	{

		// grad_F_x1 = grad_F * Grad_F_x1
		// - F1_1 * (DsInv1_1 + DsInv2_1 + DsInv3_1) - F1_2 * (DsInv1_2 + DsInv2_2 + DsInv3_2) - F1_3 * (DsInv1_3 + DsInv2_3 + DsInv3_3)
		// - F2_1 * (DsInv1_1 + DsInv2_1 + DsInv3_1) - F2_2 * (DsInv1_2 + DsInv2_2 + DsInv3_2) - F2_3 * (DsInv1_3 + DsInv2_3 + DsInv3_3)
		// - F3_1 * (DsInv1_1 + DsInv2_1 + DsInv3_1) - F3_2 * (DsInv1_2 + DsInv2_2 + DsInv3_2) - F3_3 * (DsInv1_3 + DsInv2_3 + DsInv3_3)
		FloatingType mulplier1 = DsInv(0, 0) + DsInv(1, 0) + DsInv(2, 0);
		FloatingType mulplier2 = DsInv(0, 1) + DsInv(1, 1) + DsInv(2, 1);
		FloatingType mulplier3 = DsInv(0, 2) + DsInv(1, 2) + DsInv(2, 2);

		p1Grad = {
			-FGrad(0) * mulplier1 - FGrad(1) * mulplier2 - FGrad(2) * mulplier3,
			-FGrad(3) * mulplier1 - FGrad(4) * mulplier2 - FGrad(5) * mulplier3,
			-FGrad(6) * mulplier1 - FGrad(7) * mulplier2 - FGrad(8) * mulplier3,
		};

		// grad_F_x2 = grad_F * Grad_F_x2
		// DsInv1_1 * F1_1 + DsInv1_2 * F1_2 + DsInv1_3 * F1_3
		// DsInv1_1 * F2_1 + DsInv1_2 * F2_2 + DsInv1_3 * F2_3
		// DsInv1_1 * F3_1 + DsInv1_2 * F3_2 + DsInv1_3 * F3_3
		p2Grad <<
			DsInv(0, 0) * FGrad(0) + DsInv(0, 1) * FGrad(1) + DsInv(0, 2) * FGrad(2),
			DsInv(0, 0)* FGrad(3) + DsInv(0, 1) * FGrad(4) + DsInv(0, 2) * FGrad(5),
			DsInv(0, 0)* FGrad(6) + DsInv(0, 1) * FGrad(7) + DsInv(0, 2) * FGrad(8)
			;

		// grad_F_x3 = grad_F * Grad_F_x3
		p3Grad <<
			DsInv(1, 0) * FGrad(0) + DsInv(1, 1) * FGrad(1) + DsInv(1, 2) * FGrad(2),
			DsInv(1, 0)* FGrad(3) + DsInv(1, 1) * FGrad(4) + DsInv(1, 2) * FGrad(5),
			DsInv(1, 0)* FGrad(6) + DsInv(1, 1) * FGrad(7) + DsInv(1, 2) * FGrad(8)
			;

		// grad_F_x3 = grad_F * Grad_F_x3
		p4Grad <<
			DsInv(2, 0) * FGrad(0) + DsInv(2, 1) * FGrad(1) + DsInv(2, 2) * FGrad(2),
			DsInv(2, 0)* FGrad(3) + DsInv(2, 1) * FGrad(4) + DsInv(2, 2) * FGrad(5),
			DsInv(2, 0)* FGrad(6) + DsInv(2, 1) * FGrad(7) + DsInv(2, 2) * FGrad(8)
			;
	}

	inline bool GAIA::ObjectParamsVBD::fromJson(nlohmann::json& objectParam)
	{
		ObjectParams::fromJson(objectParam);
		EXTRACT_FROM_JSON(objectParam, frictionDynamic);
		EXTRACT_FROM_JSON(objectParam, frictionEpsV);
		EXTRACT_FROM_JSON(objectParam, exponentialVelDamping);
		EXTRACT_FROM_JSON(objectParam, constantVelDamping);

		return true;
	}

	inline bool GAIA::ObjectParamsVBD::toJson(nlohmann::json& objectParam)
	{
		ObjectParams::toJson(objectParam);
		PUT_TO_JSON(objectParam, frictionDynamic);
		PUT_TO_JSON(objectParam, frictionEpsV);
		PUT_TO_JSON(objectParam, exponentialVelDamping);
		PUT_TO_JSON(objectParam, constantVelDamping);

		return true;
	}

	inline void VBDBaseTetMeshShared::initializeGPUTetMesh(VBDBaseTetMesh* pTetMesh, VBDBaseTetMeshGPU* pTetMeshGPU, bool withAccelerator)
	{
		TetMeshFEMShared::initializeGPUTetMesh(pTetMesh, pTetMeshGPU);
		activeCollisionMaskBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->activeCollisionMask.size(),
			true, pTetMesh->activeCollisionMask.data());
		activeCollisionMaskBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, activeCollisionMask);

		tetForceAndHessiansBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(12 * pTetMesh->numTets(),
			false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetForceAndHessians);
		inertiaBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->inertia.size(), true, pTetMesh->inertia.data());
		inertiaBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, inertia);
#ifdef GPU_JACOBI_DX
		positionsNewBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->inertia.size(), false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, positionsNew);
#endif // GPU_JACOBI_DX

		if (withAccelerator)
		{
			positionsPrevIterBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->inertia.size(), true);
			SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, positionsPrevIter);

		}


		collisionDataBuffer = std::make_shared<ManagedBuffer<CollisionDataGPU>>(pTetMesh->numVertices(), true);;
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, collisionData);

		pTetMeshGPU->maxVelocityMagnitude = pTetMesh->pObjParamsVBD->maxVelocityMagnitude;
		pTetMeshGPU->exponentialVelDamping = pTetMesh->pObjParamsVBD->exponentialVelDamping;
		pTetMeshGPU->constantVelDamping = pTetMesh->pObjParamsVBD->constantVelDamping;
		pTetMeshGPU->frictionDynamic = pTetMesh->pObjParamsVBD->frictionDynamic;
		pTetMeshGPU->frictionEpsV = pTetMesh->pObjParamsVBD->frictionEpsV;

	}

	inline void VBDBaseTetMeshShared::initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU)
	{
		TetMeshFEMShared::initializeGPUTetMeshForCPUDebug(pTetMesh, pTetMeshGPU);

		VBDBaseTetMeshGPU* pTetMeshGPUVBD = (VBDBaseTetMeshGPU*)pTetMeshGPU;
		VBDBaseTetMesh* pTetMeshVBD = (VBDBaseTetMesh*)pTetMesh;
		pTetMeshGPUVBD->activeCollisionMask = pTetMeshVBD->activeCollisionMask.data();
		pTetMeshGPUVBD->tetForceAndHessians = nullptr;
		pTetMeshGPUVBD->inertia = pTetMeshVBD->inertia.data();
		pTetMeshGPUVBD->collisionData = collisionDataBuffer->getCPUBuffer();

		pTetMeshGPUVBD->maxVelocityMagnitude = pTetMeshVBD->pObjParamsVBD->maxVelocityMagnitude;
		pTetMeshGPUVBD->exponentialVelDamping = pTetMeshVBD->pObjParamsVBD->exponentialVelDamping;
		pTetMeshGPUVBD->constantVelDamping = pTetMeshVBD->pObjParamsVBD->constantVelDamping;
		pTetMeshGPUVBD->frictionDynamic = pTetMeshVBD->pObjParamsVBD->frictionDynamic;
		pTetMeshGPUVBD->frictionEpsV = pTetMeshVBD->pObjParamsVBD->frictionEpsV;
	}



	inline void VBDBaseTetMeshShared::syncToGPU(bool sync, cudaStream_t stream)
	{
		TetMeshFEMShared::syncToGPU(false, stream);
		//activeCollisionMaskBuffer->toGPU(false, stream);
		collisionDataBuffer->toGPU();
		inertiaBuffer->toGPU(false, stream);

		if (sync) {
			cudaStreamSynchronize(stream);
		}
	}

	inline void VBDBaseTetMeshShared::syncToCPU(bool sync, cudaStream_t stream)
	{
		TetMeshFEMShared::syncToCPU(false, stream);
		//activeCollisionMaskBuffer->toCPU(false, stream);
		//inertiaBuffer->toCPU(false, stream);
		if (sync) {
			cudaStreamSynchronize(stream);
		}
	}

}