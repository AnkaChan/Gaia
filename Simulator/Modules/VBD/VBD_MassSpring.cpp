//#define OUTPUT_PER_ITERATION_RESULT
//#define OUTPUT_PER_ITERATION_STATISTICS
//
//#ifdef OUTPUT_PER_ITERATION_RESULT
//#include "../Framework/BasePhysicsFramework.h"
//#endif
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//#include "../Framework/BasePhysicsFramework.h"
//#include "../Timer/Timer.h"
//#endif
//
#include "VBD_MassSpring.h"
//
//#include "../Parallelization/CPUParallelization.h"
//#include <iostream>
//
//#include <MeshFrame/Utility/Str.h>
//#include <MeshFrame/Utility/IO.h>
//
//#include <Eigen/IterativeLinearSolvers>
//
//
//using namespace GAIA;
//
//void GAIA::EBDTetMesh_MassSpring::solve()
//{
//	if (pObjectParamsMaterial->localDescent)
//	{
//		solve_vertexBlockDescent();
//
//	}
//	else
//	{
//		solve_implicitEuler();
//
//	}
//
//}
//
//void GAIA::EBDTetMesh_MassSpring::solve_vertexBlockDescent()
//{
//	prevK = K;
//	prevU = U;
//
//	forwardStepInertiaAndConstantForce();
//
//	// positionsFlatten = positionsPrevFlatten;
//
//	FloatingType meritEnergy_inertia, meritEnergy_elastic;
//	FloatingType meritEnergy;
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//		std::cout << "Initial overall merit energy:" << meritEnergy << " | meritEnergy_inertia: " << meritEnergy_inertia << " | meritEnergy_elastic: " << meritEnergy_elastic << "\n";
//		}
//	);
//#if defined OUTPUT_PER_ITERATION_RESULT || defined OUTPUT_PER_ITERATION_STATISTICS
//	std::string outDebugFolder = pPhysicsFramework->outputFolder + "/Debug";
//	std::string outFolder = outDebugFolder + "/Frame_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->frameId), 4, '0');
//	MF::IO::createFolder(outFolder);
//#endif
//
//
//#ifdef OUTPUT_PER_ITERATION_RESULT
//	std::string outName = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_iteration_000.ply";
//	save(outName);
//#endif // OUTPUT_PER_ITERATION_RESULT
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//	meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//
//	std::vector<FloatingType> meritEnergy_all, meritEnergy_inertia_all, meritEnergy_elastic_all, iterationTimePerIter_all;
//
//	meritEnergy_all.push_back(meritEnergy);
//	meritEnergy_inertia_all.push_back(meritEnergy_inertia);
//	meritEnergy_elastic_all.push_back(meritEnergy_elastic);
//	double iterationTime;
//	TICK(iterationTime);
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//	FloatingType iterationTimePerIter;
//
//	for (size_t iGlobal = 0; iGlobal < pPhysicsParams->localGlobalDescentMaxIterations; iGlobal++)
//	{
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//		TICK(iterationTimePerIter);
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//
//		for (size_t iColor = 0; iColor < verticesColoring.size(); iColor++)
//		{
//			const std::vector<int> &color = verticesColoring[iColor];
//			cpu_parallel_for(0, color.size(), [&](int iV) {
//				localMeritEnergyVertexBlockCoordinateDesent(color[iV]);
//			});
//		}
//		//for (size_t iV = 0; iV < numVertices(); iV++)
//		//{
//		//	localMeritEnergyVertexBlockCoordinateDesent(iV);
//		//}
//		debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//			meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//			std::cout << "Merit energy after iter " << iGlobal << " :" << meritEnergy << " | meritEnergy_inertia: " << meritEnergy_inertia << " | meritEnergy_elastic: " << meritEnergy_elastic << "\n";
//			});
//		//energyProjectionGlobal();
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//		TOCK(iterationTimePerIter);
//		iterationTimePerIter_all.push_back(iterationTimePerIter);
//
//		meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//		meritEnergy_all.push_back(meritEnergy);
//		meritEnergy_inertia_all.push_back(meritEnergy_inertia);
//		meritEnergy_elastic_all.push_back(meritEnergy_elastic);
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//
//#ifdef OUTPUT_PER_ITERATION_RESULT
//		outName = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_iteration_" + MF::STR::padToNew(std::to_string(iGlobal + 1), 3, '0') + ".ply";
//		save(outName);
//#endif
//	}
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//	TOCK(iterationTime);
//	nlohmann::json outConvergenceStats;
//	outConvergenceStats["iterations"] = pPhysicsParams->localGlobalDescentMaxIterations;
//	outConvergenceStats["time"] = iterationTime;
//	outConvergenceStats["timePerIter"] = iterationTimePerIter_all;
//	outConvergenceStats["meritEnergy"] = meritEnergy_all;
//	outConvergenceStats["meritEnergy_inertia"] = meritEnergy_inertia_all;
//	outConvergenceStats["meritEnergy_elastic"] = meritEnergy_elastic_all;
//	std::string outStatsFile = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_Stats.json" ;
//	std::ofstream ofs(outStatsFile);
//	if (ofs.is_open())
//	{
//		ofs << outConvergenceStats.dump(2);
//	}
//	ofs.close();
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//	/*energyResidual = K + U - prevK - prevU;
//	std::cout << "Energy Residual after local steps: " 
//		<< K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";*/
//
//	// energyConservationLineSearch_dualSpaceStepping(energyResidual, );
//
//	// energyConjugateGradientDescent_dualSpace();
//
//	// energy projection
//	//energyConservationLineSearch_dualSpaceStepping();
//
//	//forwardStepSymplecticEuler();
//	//if (pPhysicsParams->applyEnergyLineSearch)
//	//{
//	//	positionsFlatten = positionsPrevFlatten;
//	//	for (size_t iV = 0; iV < numVertices(); iV++)
//	//	{
//	//		if (!fixedMask(iV))
//	//		{
//	//			energyConservationLineSearch_dualSpaceStepping_perVertex(iV);
//	//		}
//	//	}
//	//}
//
//	// Evaluate per vertex and per element energy
//
//	U = evaluatePotentialEnergy(positionsFlatten, U_perElement);
//	K = evaluateKinecticEnergy(velocityFlatten, K_perVertex);
//
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		K = evaluateKinecticEnergy(velocityFlatten);
//		FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//		U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//
//		FloatingType energyResidual = K + U - prevK - prevU;
//		std::cout << "Final Energy Residual: " << K + U - prevK - prevU << " | U: " << U << "| K: " << K <<  "\n";
//		});
//
//}
//
//void GAIA::EBDTetMesh_MassSpring::solve_implicitEuler()
//{
//	prevK = K;
//	prevU = U;
//
//	//hessianAssembler.initializeForMassSpring(numVertices(), edges);
//
//	if (hessianAssembler.hessian.rows() == 0)
//	{
//		hessianAssembler.initializeForMassSpring(numVertices(), edges, pObjectParamsMaterial->fixedPoints);
//
//		hessianAssembler.toEigenTriples(true);
//		hessianAssembler.toSparseMat();
//	}
//
//	forwardStepInertiaAndConstantForce();
//
//	// positionsFlatten = positionsPrevFlatten;
//
//	FloatingType meritEnergy_inertia, meritEnergy_elastic;
//	FloatingType meritEnergy;
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//		std::cout << "Initial overall merit energy:" << meritEnergy << " | meritEnergy_inertia: " << meritEnergy_inertia << " | meritEnergy_elastic: " << meritEnergy_elastic << "\n";
//		}
//	);
//#if defined OUTPUT_PER_ITERATION_RESULT || defined OUTPUT_PER_ITERATION_STATISTICS
//	std::string outDebugFolder = pPhysicsFramework->outputFolder + "/Debug";
//	std::string outFolder = outDebugFolder + "/Frame_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->frameId), 4, '0');
//	MF::IO::createFolder(outFolder);
//#endif
//
//
//#ifdef OUTPUT_PER_ITERATION_RESULT
//	std::string outName = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_iteration_000.ply";
//	save(outName);
//#endif // OUTPUT_PER_ITERATION_RESULT
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//	meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//
//	std::vector<FloatingType> meritEnergy_all, meritEnergy_inertia_all, meritEnergy_elastic_all, iterationTimePerIter_all;
//
//	meritEnergy_all.push_back(meritEnergy);
//	meritEnergy_inertia_all.push_back(meritEnergy_inertia);
//	meritEnergy_elastic_all.push_back(meritEnergy_elastic);
//	double iterationTime;
//	TICK(iterationTime);
//	FloatingType iterationTimePerIter;
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//
//	for (size_t iGlobal = 0; iGlobal < pPhysicsParams->localGlobalDescentMaxIterations; iGlobal++)
//	{
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//		TICK(iterationTimePerIter);
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//
//		solve_implicitEuler_iteration();
//
//		debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//			meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//			std::cout << "Merit energy after iter " << iGlobal << " :" << meritEnergy << " | meritEnergy_inertia: " << meritEnergy_inertia << " | meritEnergy_elastic: " << meritEnergy_elastic << "\n";
//			});
//		//energyProjectionGlobal();
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//		TOCK(iterationTimePerIter);
//		iterationTimePerIter_all.push_back(iterationTimePerIter);
//
//		meritEnergy = evaluateMeritEnergy(meritEnergy_inertia, meritEnergy_elastic);
//		meritEnergy_all.push_back(meritEnergy);
//		meritEnergy_inertia_all.push_back(meritEnergy_inertia);
//		meritEnergy_elastic_all.push_back(meritEnergy_elastic);
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//
//#ifdef OUTPUT_PER_ITERATION_RESULT
//		outName = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_iteration_" + MF::STR::padToNew(std::to_string(iGlobal + 1), 3, '0') + ".ply";
//		save(outName);
//#endif
//	}
//
//	updateVelocityBasedOnPositionChange();
//
//
//#ifdef OUTPUT_PER_ITERATION_STATISTICS
//	TOCK(iterationTime);
//	nlohmann::json outConvergenceStats;
//	outConvergenceStats["iterations"] = pPhysicsParams->localGlobalDescentMaxIterations;
//	outConvergenceStats["time"] = iterationTime;
//	outConvergenceStats["timePerIter"] = iterationTimePerIter_all;
//	outConvergenceStats["meritEnergy"] = meritEnergy_all;
//	outConvergenceStats["meritEnergy_inertia"] = meritEnergy_inertia_all;
//	outConvergenceStats["meritEnergy_elastic"] = meritEnergy_elastic_all;
//	std::string outStatsFile = outFolder + "/" + +"Step_" + MF::STR::padToNew(std::to_string(pPhysicsFramework->substep), 4, '0') + "_Stats.json";
//	std::ofstream ofs(outStatsFile);
//	if (ofs.is_open())
//	{
//		ofs << outConvergenceStats.dump(2);
//	}
//	ofs.close();
//#endif // OUTPUT_PER_ITERATION_STATISTICS
//	/*energyResidual = K + U - prevK - prevU;
//	std::cout << "Energy Residual after local steps: "
//		<< K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";*/
//
//		// energyConservationLineSearch_dualSpaceStepping(energyResidual, );
//
//		// energyConjugateGradientDescent_dualSpace();
//
//		// energy projection
//		//energyConservationLineSearch_dualSpaceStepping();
//
//		//forwardStepSymplecticEuler();
//		//if (pPhysicsParams->applyEnergyLineSearch)
//		//{
//		//	positionsFlatten = positionsPrevFlatten;
//		//	for (size_t iV = 0; iV < numVertices(); iV++)
//		//	{
//		//		if (!fixedMask(iV))
//		//		{
//		//			energyConservationLineSearch_dualSpaceStepping_perVertex(iV);
//		//		}
//		//	}
//		//}
//
//		// Evaluate per vertex and per element energy
//
//	U = evaluatePotentialEnergy(positionsFlatten, U_perElement);
//	K = evaluateKinecticEnergy(velocityFlatten, K_perVertex);
//
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		K = evaluateKinecticEnergy(velocityFlatten);
//		FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//		U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//
//		FloatingType energyResidual = K + U - prevK - prevU;
//		std::cout << "Final Energy Residual: " << K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";
//		});
//
//}
//
//void GAIA::EBDTetMesh_MassSpring::solve_implicitEuler_iteration()
//{
//	FloatingType dtSqrReciprocal = 1.f / (pPhysicsParams->dt * pPhysicsParams->dt);
//	
//	for (int k = 0; k < hessianAssembler.hessian.outerSize(); ++k) {
//		for (Eigen::SparseMatrix<FloatingType>::InnerIterator it(hessianAssembler.hessian, k); it; ++it)
//		{
//			it.valueRef() = 0.f;
//		}
//	}
//	internalForce.setZero();
//
//	VecDynamic gradient;
//	gradient.resize(hessianAssembler.numVerts * 3);
//	gradient.setZero();
//
//	FloatingType timeConsumptionAssemble, timeConsumptionSolve;
//
//	TICK(timeConsumptionAssemble);
//	for (size_t v1 = 0; v1 < numVertices(); v1++)
//	{
//		if (fixedMask[v1]) {
//			continue;
//		}
//		int hesRefiV1 = hessianAssembler.getHessianRef(v1);
//
//		// evaluate inertia energy's Hessian
//		hessianAssembler.hessian.coeffRef(hesRefiV1,   hesRefiV1)   = vertexMass(v1) * dtSqrReciprocal;
//		hessianAssembler.hessian.coeffRef(hesRefiV1+1, hesRefiV1+1) = vertexMass(v1) * dtSqrReciprocal;
//		hessianAssembler.hessian.coeffRef(hesRefiV1+2, hesRefiV1+2) = vertexMass(v1) * dtSqrReciprocal;
//
//		Vec3 VInertia = BLOCK_DIM3x1_AT(inertia, v1);
//
//		BLOCK_DIM3x1_AT(gradient, hessianAssembler.vertexIdConverter[v1]) = (vertex(v1) - VInertia) * vertexMass(v1) * dtSqrReciprocal;
//	}
//	
//	// evaluate elastic energy's Hessian and force
//	for (size_t iSpring = 0; iSpring < edges.cols(); iSpring++)
//	{
//		IdType v1 = edges(0, iSpring);
//		IdType v2 = edges(1, iSpring);
//		Vec3 diff = mVertPos.col(v1) - mVertPos.col(v2);
//		FloatingType l = diff.norm();
//
//		FloatingType l0 = orgLengths(iSpring);
//
//		Vec3 force = (adjustedSpringStiffness(iSpring) * (l0 - l) / l) * diff;
//
//		// evaluate hessian
//		Mat3 h_1_1;
//		h_1_1 = adjustedSpringStiffness(iSpring) * (Mat3::Identity() - (l0 / l) * (Mat3::Identity() - diff * diff.transpose() / (l * l)));
//
//		// std::cout << "Hessian for edge: (" << v1 << ", " << v2 << ") is:\n"
//		// 	<< h_1_1 << "\n";
//		//std::cout << "Hessian (numerical) for edge: (" << v1 << ", " << v2 << ") is:\n"
//
//		//	<< computeNumericalHessianForSpring(iSpring) << "\n";
//
//		int hesRefV1 = hessianAssembler.getHessianRef(v1);
//		int v1Filtered = hessianAssembler.vertexIdConverter[v1];
//		int hesRefV2 = hessianAssembler.getHessianRef(v2);
//		int v2Filtered = hessianAssembler.vertexIdConverter[v2];
//		Mat3 h_1_1_n = -h_1_1;
//
//		if (!fixedMask(v1))
//		{
//			BLOCK_DIM3x1_AT(gradient, v1Filtered) -= force;
//			
//			hessianAssembler.addToHessianBlock(hesRefV1, hesRefV1, h_1_1);
//			if (v2Filtered >= 0)
//			{
//				hessianAssembler.addToHessianBlock(hesRefV1, hesRefV2, h_1_1_n);
//			}
//		}
//
//		if (!fixedMask(v2)) 
//		{
//			BLOCK_DIM3x1_AT(gradient, v2Filtered) += force;
//			hessianAssembler.addToHessianBlock(hesRefV2, hesRefV2, h_1_1);
//			if (v1Filtered >= 0)
//			{
//				hessianAssembler.addToHessianBlock(hesRefV2, hesRefV1, h_1_1_n);
//			}
//		}
//	}
//	TOCK(timeConsumptionAssemble);
//
//	TICK(timeConsumptionSolve);
//	//Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<FloatingType>> solver;
//	Eigen::ConjugateGradient<Eigen::SparseMatrix<FloatingType>> solver;
//	solver.setMaxIterations(200);
//	solver.setTolerance(1e-5f);
//	solver.compute(hessianAssembler.hessian);
//	//std::cout << "internalForce: " << internalForce.transpose() << std::endl;
//	//std::cout << "gradientInertia: " << gradientInertia.transpose() << std::endl;
//	//std::cout << "hessian: \n" << Eigen::MatrixXf(hessianAssembler.hessian) << std::endl;
//	//Eigen::MatrixXf hDense(hessianAssembler.hessian);
//	//Eigen::BDCSVD<Eigen::MatrixXf> svd(hDense, Eigen::ComputeThinU | Eigen::ComputeThinV);
//	//VecDynamic eigenVals = svd.singularValues();;
//	//std::cout << "Eigen values of Hessian:\n" << eigenVals.transpose() << std::endl;
//
//	//Eigen::SparseLU<Eigen::SparseMatrix<FloatingType>, Eigen::COLAMDOrdering<int>>  solver;
//	// fill A and b;
//	// Compute the ordering permutation vector from the structural pattern of A
//	//solver.analyzePattern(hessianAssembler.hessian);
//	// Compute the numerical factorization 
//	//solver.factorize(hessianAssembler.hessian);
//	//Use the factors to solve the linear system 
//
//	VecDynamic dx = solver.solve(gradient);
//
//	VecDynamic dxOrg = hessianAssembler.getOriginalVector(dx);
//	TOCK(timeConsumptionSolve);
//
//	std::cout << "#iterations:     " << solver.iterations() << std::endl;
//	std::cout << "estimated error: " << solver.error() << std::endl;
//	std::cout << "timeConsumptionAssemble: " << timeConsumptionAssemble << " | " << timeConsumptionSolve << std::endl;
//
//	positionsFlatten -= dxOrg;
//
//}
//
//void GAIA::EBDTetMesh_MassSpring::solve_localEnergyConsistency()
//{
//	prevK = K;
//	prevU = U;
//
//	prevK_perVertex = K_perVertex;
//	prevU_perElement = U_perElement;
//
//	forwardStepSymplecticEuler();
//	K = evaluateKinecticEnergy(velocityFlatten);
//	FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//	U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//	FloatingType energyResidual = K + U - prevK - prevU;
//	std::cout << "Initial Energy Residual LS: " << K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";
//
//	if (energyResidual > pPhysicsParams->globalEnergyProjectionConvergeThres)
//	{
//		energyProjectionLocal();
//	}
//
//	energyProjectionGlobal();
//
//	/*energyResidual = K + U - prevK - prevU;
//	std::cout << "Energy Residual after local steps: "
//		<< K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";*/
//
//		// energyConservationLineSearch_dualSpaceStepping(energyResidual, );
//
//		// energyConjugateGradientDescent_dualSpace();
//
//		// energy projection
//		//energyConservationLineSearch_dualSpaceStepping();
//
//		//forwardStepSymplecticEuler();
//		//if (pPhysicsParams->applyEnergyLineSearch)
//		//{
//		//	positionsFlatten = positionsPrevFlatten;
//		//	for (size_t iV = 0; iV < numVertices(); iV++)
//		//	{
//		//		if (!fixedMask(iV))
//		//		{
//		//			energyConservationLineSearch_dualSpaceStepping_perVertex(iV);
//		//		}
//		//	}
//		//}
//
//		// Evaluate per vertex and per element energy
//
//	U = evaluatePotentialEnergy(positionsFlatten, U_perElement);
//	K = evaluateKinecticEnergy(velocityFlatten, K_perVertex);
//
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		K = evaluateKinecticEnergy(velocityFlatten);
//		FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//		U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//
//		FloatingType energyResidual = K + U - prevK - prevU;
//		std::cout << "Final Energy Residual: " << K + U - prevK - prevU << " | U: " << U << "| K: " << K << "\n";
//		});
//}
//
//
//void GAIA::EBDTetMesh_MassSpring::solve_energyProjection()
//{
//	// FloatingType prevE = evaluateInternalPotentialEnergy(positionsFlatten);
//
//	//FloatingType dt = 1e-5f;
//	//VecDynamic internalForceFlattenNew;
//	//VecDynamic internalForceFlattenPrev_noGrav;
//	//internalForceFlattenNew.resizeLike(internalForceFlatten);
//	//internalForceFlattenPrev_noGrav.resizeLike(internalForceFlatten);
//	//evaluateInternalForce(positionsFlatten, internalForceFlattenPrev_noGrav,false);
//	//for (size_t i = 0; i < 3*numVertices(); i++)
//	//{
//	//	VecDynamic newPosFlatten = 1.f*positionsFlatten;
//
//	//	newPosFlatten(i) += dt;
//
//	//	FloatingType newE = evaluateInternalPotentialEnergy(newPosFlatten);
//	//	evaluateInternalForce(newPosFlatten, internalForceFlattenNew, false);
//	//	//std::cout << "Force numerical: " << -(newE - prevE) / dt  << " | "
//	//	//	<< "Force analytical: " << internalForceFlattenNew(i) << "\n";
//	//	assert(abs((newE - prevE)/dt + (internalForceFlattenPrev_noGrav(i) + internalForceFlattenNew(i)) * 0.5 ) < 1e-4);
//	//}
//	evaluateInternalForce(positionsFlatten, internalForceFlatten);
//	handleForceConstraints();
//
//	prevK = K;
//	prevU = U;
//
//	//forwardStepSymplecticEuler();
//	forwardStepForwardEuler();
//
//	// evaluateInternalForce();
//	// evaluate energy residual
//
//
//	K = evaluateKinecticEnergy(velocityFlatten);
//	FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//	U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//
//	FloatingType energyResidual = K + U - prevK - prevU;
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		std::cout << "Initial Energy Residual: " << energyResidual << " | Graviational energy: " << graviationalEnergy
//			<< " | Elasitic energy: " << U - graviationalEnergy << " | Kinietic energy: " << K << "\n";
//		});
//
//	//
//	////// energy projection
//	// energyConservationLineSearchDiagonalPreconditioner(energyResidual, prevK + prevU);
//	//energyConservationPotentialOnlyLineSearch_KBounded(energyResidual, prevK + prevU);
//	energyConservationPotentialOnlyLineSearch_KUBalanced(energyResidual, prevK + prevU);
//
//	K = evaluateKinecticEnergy(velocityFlatten);
//	U = evaluateInternalPotentialEnergy(positionsFlatten) + evaluateGravitationalEnergy(positionsFlatten);
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		std::cout << "Final Energy Residual: " << K + U - prevK - prevU << "\n";
//		});
//}
//
//void GAIA::EBDTetMesh_MassSpring::energyProjectionLocal()
//{
//	for (size_t iGlobal = 0; iGlobal < pPhysicsParams->localEnergyProjectionMaxIterOutterLoop; iGlobal++)
//	{
//		for (size_t iV = 0; iV < numVertices(); iV++)
//		{
//			energyProjectionLocalForVertex(iV);
//		}
//	}
//}
//
//void GAIA::EBDTetMesh_MassSpring::energyProjectionGlobal()
//{
//
//	FloatingType lambda_l = 0, res_l = -prevK;
//
//	/*if (!pPhysicsParams->applyEnergyLineSearch)
//	{
//		velocityFlatten = velocityFlatten * sqrt(pObjParamsVBD->exponentialVelDamping);
//		return;
//	}*/
//
//	K = evaluateKinecticEnergy(velocityFlatten);
//	FloatingType graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//	U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//	// evaluateVelocityNorms(velocityFlatten, velNorms);
//
//	// Utility::printVecInfo(velocityFlatten, "velocityFlatten");
//
//	FloatingType KOrg = K;
//
//	FloatingType lambda_r = pPhysicsParams->dt, res_r = K + U - prevK - prevU;
//	FloatingType d_in = -1;
//	FloatingType lambda = lambda_r;
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		std::cout << "Initial Energy Residual GS: " << res_r << " | Graviational energy: " << graviationalEnergy
//			<< " | Elasitic energy: " << U - graviationalEnergy << " | Kinietic energy: " << K << "\n";
//		});
//
//	if (res_r < pPhysicsParams->linesearchConvergeThres) {
//		// TODO: enlarge lambda_r or store energy debt
//		res_l = res_r;
//		lambda_l = lambda_r;
//		
//		while (res_r < 0)
//		{
//			lambda_r *= pPhysicsParams->linesearchIntervalGrowthMul;
//			K = pow(lambda_r / pPhysicsParams->dt, 2) * KOrg;
//			positionsFlatten = positionsPrevFlatten + lambda_r * velocityFlatten;
//
//			graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//			U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//			res_r = K + U - prevU - prevK;
//		}
//		lambda = lambda_r;
//	}
//	FloatingType res = res_r;
//
//	int iter = 0;
//
//	while (abs(res) > pPhysicsParams->linesearchConvergeThres
//		&& lambda_r - lambda_l > pPhysicsParams->linesearchIntervalSizeThres
//		) {
//		evaluateInternalForce(positionsFlatten, internalForceFlatten, true);
//		handleForceConstraints();
//		getPositionsBlock(gradOrg) = -internalForceFlatten - externalForceFlatten;
//
//		FloatingType dL_dLambda = getPositionsBlock(gradOrg).dot(velocityFlatten) + 2*lambda * KOrg;
//		FloatingType d_lambda = -res / dL_dLambda;
//		lambda += d_lambda;
//
//		if (abs(d_lambda) < pPhysicsParams->linesearchIntervalSizeThres 
//			|| lambda < lambda_l || lambda > lambda_r) {
//			// desecent direction points out of the interval, doing linear intepolation
//			FloatingType t = res_r / (res_r - res_l);
//
//			lambda = lambda_l * t + lambda_r * (1 - t);
//		}
//
//		positionsFlatten = positionsPrevFlatten + lambda * velocityFlatten;
//
//		K = pow(lambda / pPhysicsParams->dt, 2) * KOrg;
//		graviationalEnergy = evaluateGravitationalEnergy(positionsFlatten);
//		U = evaluateInternalPotentialEnergy(positionsFlatten) + graviationalEnergy;
//
//		res = K + U - prevU - prevK;
//
//		if (res < 0)
//		{
//			lambda_l = lambda;
//			res_l = res;
//		}
//		else {
//			lambda_r = lambda;
//			res_r = res;
//		}
//		debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG_VEBOSE, [&]() {
//			std::cout << "GS Iteration: " << iter << " energy residual: " << res << " lambda: " << lambda
//				<< " d_lambda: " << d_lambda << " l_l: " << lambda_l << " l_l: " << lambda_r
//				//<< "\nGraviational energy: " << graviationalEnergy
//				//<< " | Elasitic energy: " << U - graviationalEnergy << " | Kinietic energy: " << K << "\n"
//				// << "position v 200:" << vertex(200).transpose() 
//				<< "\n"
//				;
//
//			});
//
//		iter++;
//	}
//	velocityFlatten = velocityFlatten * (lambda / pPhysicsParams->dt);
//	// velocityFlatten = velocityFlatten * (lambda / pPhysicsParams->dt) * sqrt(pObjParamsVBD->exponentialVelDamping);
//	// K = K * pObjParamsVBD->exponentialVelDamping;
//
//	debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG, [&]() {
//		std::cout << "Final Energy Residual GS: " << res << " | Graviational energy: " << graviationalEnergy
//			<< " | Elasitic energy: " << U - graviationalEnergy << " | Kinietic energy: " << K << "\n";
//		});
//}
//
//void GAIA::EBDTetMesh_MassSpring::energyProjectionLocalForVertex(int iV)
//{
//	if (fixedMask[iV])
//	{
//		return;
//	}
//
//	FloatingType U_prev_div = 0.f;
//	FloatingType U_div = 0.f;
//	FloatingType K_prev = prevK_perVertex[iV];
//	FloatingType K = evaluateKinaticEnergyPerVertex(iV);
//	FloatingType K2 = evaluateKinaticEnergyPerVertexFromPosDiff(iV);
//
//	for (size_t iNeiSpring = 0; iNeiSpring < vertexNeighborEdges[iV].size(); iNeiSpring++)
//	{
//		IdType neiSpringId = vertexNeighborEdges[iV][iNeiSpring];
//
//		U_prev_div += prevU_perElement[neiSpringId];
//	}
//
//	U_prev_div = U_prev_div * massDivisionRatio - vertexPrevPos(iV)(GRAVITY_AXIS) * pPhysicsParams->gravity(GRAVITY_AXIS) * vertexMass(iV);
//	U_div = evaluateLocalInternalEnergyPerVertex(iV) * massDivisionRatio + evaluateGraviationalEnergyPerVertex(iV);
//
//	FloatingType prevEnergy = K_prev + U_prev_div;
//
//	// res = U_div - U_prev_div + (K - K_prev)
//	FloatingType res = U_div + K - prevEnergy;
//	int iter = 0;
//	Vec3 v_force;
//	Vec3 U_grad;
//	Vec3 K_grad;
//	Vec3 grad;
//	Vec3 bestPos = vertex(iV);
//	FloatingType bestRes = res;
//	while (iter < pPhysicsParams->localEnergyProjectionMaxIterInnerLoop
//		&& abs(res) > pPhysicsParams->localEnergyProjectionConvergeThres)
//	{
//		if (iter == 0)
//		{
//			//std::cout << "Local energy residual: " << res << "\n";
//		}
//
//		evaluateInternalForcePerVertex(iV, v_force, false);
//		U_grad = -massDivisionRatio * v_force;
//		Vec3 gravity = -pPhysicsParams->gravity * vertexMass(iV);
//		U_grad += gravity;
//		K_grad = (vertex(iV) - vertexPrevPos(iV)) * vertexMass(iV) / (pPhysicsParams->dt * pPhysicsParams->dt);
//
//		grad = U_grad * (U_div-U_prev_div)/(res) + K_grad * (K-K_prev) / (res);
//		// grad = U_grad;
//
//
//		// newton step
//		FloatingType lambda = -1.f * res / grad.squaredNorm();
//
//		vertex(iV) += grad * lambda;
//
//		// re-evaluate grad
//		U_div = evaluateLocalInternalEnergyPerVertex(iV) * massDivisionRatio + evaluateGraviationalEnergyPerVertex(iV);
//		K = evaluateKinaticEnergyPerVertexFromPosDiff(iV);
//
//		res = U_div + K - prevEnergy;
//		// std::cout << "Local energy residual: " << res << ", best: " << bestRes << "\n";
//
//		if (abs(res) < bestRes)
//		{
//			bestPos = vertex(iV);
//			bestRes = abs(res);
//		}
//
//		iter++;
//	}
//
//	vertex(iV) = bestPos;
//
//	velocity(iV) = (bestPos - vertexPrevPos(iV)) / (pPhysicsParams->dt);
//}
//
//void GAIA::EBDTetMesh_MassSpring::localMeritEnergyVertexBlockCoordinateDesent(int iV)
//{
//	if (fixedMask[iV])
//	{
//		return;
//	}
//
//	Vec3 VInertia = BLOCK_DIM3x1_AT(inertia, iV);
//
//	FloatingType dtSqrReciprocal = 1.f / (pPhysicsParams->dt * pPhysicsParams->dt);
//
//	// evaluate hessian
//	Mat3 h;
//	Vec3 force;
//	evaluateInternalForcePerVertex(iV, force, false);
//	// inertiaEnergy = vertexMass(iV) * (vertex(iV) - BLOCK_DIM3x1_AT(inertia, iV)).squaredNorm() / (dt2SqrReciprocal);
//	Vec3 grad = vertexMass(iV) * (vertex(iV) - VInertia) * (dtSqrReciprocal) - force;
//
//	FloatingType e0 = evaluateMeritEnergyPervertex(iV);
//
//	FloatingType stepLength = 1.0;
//	FloatingType stepLengthReduceRate = 1.0f;
//
//	int iIter = 0;
//	while (grad.squaredNorm() > pow(pPhysicsParams->localMeritEnergyDescentConvergeThres, 2)
//		//&& e0 > pPhysicsParams->localMeritEnergyDescentConvergeThres
//		&& iIter < pPhysicsParams->localMeritEnergyDescentMaxIterOutterLoop
//		)
//	{
//		debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG_VEBOSE_2, [&]() {
//			std::cout << "Initial merit energy for Vertex " << iV << ":" << e0 << " | norm of gradient: " << grad.norm() << "\n";
//			}
//		);
//
//		evaluateLocalInternalEnergyHessianPerVertex(iV, h);
//		//h = Mat3::Zero();
//
//		if (pPhysicsParams->filterHessian)
//		{
//			Eigen::JacobiSVD<Mat3> svd(h, Eigen::ComputeFullU | Eigen::ComputeFullV);
//
//			Vec3 eigenVals = svd.singularValues();
//
//			Mat3 diagonalEV = Mat3::Zero();
//			for (size_t iDim = 0; iDim < 3; iDim++)
//			{
//				if (eigenVals(iDim) > 0) diagonalEV(iDim, iDim) = eigenVals(iDim);
//				else
//				{
//					std::cout << "Negative Hessian ecountered: " << eigenVals.transpose() << "\n";
//				}
//			}
//
//			h = svd.matrixU() * diagonalEV * svd.matrixV().transpose();
//		}
//
//		h += vertexMass(iV) * dtSqrReciprocal * Mat3::Identity();
//		Vec3 d_x = -h.colPivHouseholderQr().solve(grad);
//
//		Vec3 vPrev = BLOCK_DIM3x1_AT(positionsFlatten, iV);
//		//std::cout << "V before linesearch: " << vPrev.transpose() << "\n";
//
//		FloatingType alpha = localMeritEnergyVertexBlockCoordinateBackTrackingLineSearch(iV, e0, d_x, stepLength, 0.0f, 0.5f, pPhysicsParams->localMeritEnergyDescentMaxIterInnerLoop);
//		iIter++;
//		stepLength *= stepLengthReduceRate;
//		//std::cout << "V afer linesearch: " << BLOCK_DIM3x1_AT(positionsFlatten, iV).transpose() << "\n";
//		//std::cout << "Change: " << (BLOCK_DIM3x1_AT(positionsFlatten, iV) - vPrev).transpose() << "\n";
//
//		evaluateInternalForcePerVertex(iV, force, false);
//
//		Vec3 diffInertia = (vertex(iV) - VInertia);
//		grad = vertexMass(iV) * diffInertia * (dtSqrReciprocal) - force;
//		e0 = evaluateMeritEnergyPervertex(iV);
//
//		debugInfoGen(pPhysicsParams->debugVerboseLvl, DEBUG_LVL_DEBUG_VEBOSE_2, [&]() {
//			std::cout << "Merit energy for Vertex " << iV << " at iter " << iIter << " :" << e0 
//				<< " | norm of gradient: " << grad.norm() 
//				<< " | norm of diffInertia: " << diffInertia.norm()
//				<< " | norm of force: " << force.norm()
//				<<"\n";
//			});
//	}
//
//	velocity(iV) = (vertex(iV) - vertexPrevPos(iV)) / pPhysicsParams->dt;
//
//	//FloatingType meritEnergy = ;
//	//while (iter < pPhysicsParams->localEnergyProjectionMaxIterInnerLoop
//	//	&& abs(res) > pPhysicsParams->localEnergyProjectionConvergeThres)
//	//{
//	//	if (iter == 0)
//	//	{
//	//		//std::cout << "Local energy residual: " << res << "\n";
//	//	}
//
//	//	evaluateInternalForcePerVertex(iV, v_force, false);
//	//	U_grad = -massDivisionRatio * v_force;
//	//	Vec3 gravity = -pPhysicsParams->gravity * vertexMass(iV);
//	//	U_grad += gravity;
//	//	K_grad = (vertex(iV) - vertexPrevPos(iV)) * vertexMass(iV) / (pPhysicsParams->dt * pPhysicsParams->dt);
//
//	//	grad = U_grad * (U_div - U_prev_div) / (res)+K_grad * (K - K_prev) / (res);
//	//	// grad = U_grad;
//
//
//	//	// newton step
//	//	FloatingType lambda = -1.f * res / grad.squaredNorm();
//
//	//	vertex(iV) += grad * lambda;
//
//	//	// re-evaluate grad
//	//	U_div = evaluateLocalInternalEnergyPerVertex(iV) * massDivisionRatio + evaluateGraviationalEnergyPerVertex(iV);
//	//	K = evaluateKinaticEnergyPerVertexFromPosDiff(iV);
//
//	//	res = U_div + K - prevEnergy;
//	//	// std::cout << "Local energy residual: " << res << ", best: " << bestRes << "\n";
//
//	//	if (abs(res) < bestRes)
//	//	{
//	//		bestPos = vertex(iV);
//	//		bestRes = abs(res);
//	//	}
//
//	//	iter++;
//	//}
//
//	//vertex(iV) = bestPos;
//
//	//velocity(iV) = (bestPos - vertexPrevPos(iV)) / (pPhysicsParams->dt);
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::localMeritEnergyVertexBlockCoordinateBackTrackingLineSearch(int iV, FloatingType E0, Vec3 dx, 
//	FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters)
//{
//	Vec3 x_0 = vertex(iV);
//
//	FloatingType bestAlpha = 0.f;
//	FloatingType bestEnergy = E0;
//
//	// dx.dot(grad), since 
//	FloatingType m = dx.squaredNorm();
//
//
//	for (size_t iIter = 0; iIter < maxNumIters; iIter++)
//	{
//		vertex(iV) = alpha * dx + x_0;
//
//		FloatingType e = evaluateMeritEnergyPervertex(iV);
//
//
//		if (e < bestEnergy)
//		{
//			bestAlpha = alpha;
//			bestEnergy = e;
//		}
//
//		// first Wolfe condition 
//		if (e < E0 + alpha * c * m)
//		{
//			break;
//		}
//		else
//		{
//			alpha = alpha * tau;
//		}
//	}
//
//	vertex(iV) = bestAlpha *dx + x_0;
//
//	return bestAlpha;
//}
//
//void GAIA::EBDTetMesh_MassSpring::initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF, 
//	VBDPhysicsParameters::SharedPtr inPhysicsParams, BasePhysicFramework* in_pPhysicsFramework)
//{
//	pObjectParamsMaterial = std::static_pointer_cast<ObjectParametersEBDMassSpring>(inMaterialParams);
//	VBDBaseTetMesh::initialize(inMaterialParams, pTM_MF, inPhysicsParams, in_pPhysicsFramework);
//
//	nEdges = pTM_MF->numEdges();
//
//	orgLengths.resize(nEdges);
//	edges.resize(2, nEdges);
//
//	vertexNeighborEdges.resize(numVertices());
//
//	int iEdge = 0;
//	for (TetMeshMF::EPtr pE : TIt::TM_EIterator(pTM_MF.get()))
//	{
//		edges.col(iEdge) << pE->vertex1()->id(), pE->vertex2()->id();
//		orgLengths(iEdge) = (mVertPos.col(pE->vertex1()->id()) - mVertPos.col(pE->vertex2()->id())).norm();
//
//		int32_t* edge = edges.data() + 2 * iEdge;
//		assert(pE->vertex1()->id() == edge[0] && pE->vertex2()->id() == edge[1]);
//
//		vertexNeighborEdges[edge[0]].push_back(iEdge);
//		vertexNeighborEdges[edge[1]].push_back(iEdge);
//
//		++iEdge;
//	}
//	U_perElement.resize(nEdges);
//	avgSpringLength = orgLengths.mean();
//
//	//externalForce = VecDynamic::Zero(numVertices() * 3);
//
//	//springEnergy.resize(nEdges);
//
//	//dv.resizeLike(internalForce);
//	//b.resizeLike(internalForce);
//	//f0.resizeLike(verlocityFlatten);
//	//invMass3x.resizeLike(internalForce);
//	//for (size_t iV = 0; iV < vertexInvMass.size(); iV++)
//	//{
//	//	invMass3x(iV * 3) = vertexInvMass(iV);
//	//	invMass3x(iV * 3 + 1) = vertexInvMass(iV);
//	//	invMass3x(iV * 3 + 2) = vertexInvMass(iV);
//	//}
//
//	if (pObjectParamsMaterial->lengthBasedMass)
//	{
//		rebalanceMassBasedOnSpringLength();
//	}
//
//	std::cout << "Max vertex mass: " << vertexMass.maxCoeff() << " | min vertex mass: " << vertexMass.minCoeff() << " | mean vertex mass: " << vertexMass.mean() << "\n";
//	int minMassVertID;
//	vertexMass.minCoeff(&minMassVertID);
//	std::cout << "Min mass vertex id: " << minMassVertID << "\n";
//
//	std::cout << "Max spring length: " << orgLengths.maxCoeff() << " | min spring length: " << orgLengths.minCoeff() << " | mean spring length: " << vertexMass.mean() << "\n";
//
//}
//
//void GAIA::EBDTetMesh_MassSpring::rebalanceMassBasedOnSpringLength()
//{
//	FloatingType totalMass = vertexMass.sum();
//	FloatingType totalSpringLen = orgLengths.sum();
//
//	VecDynamic vertexAdjacentELengthsAll = decltype(vertexMass)::Zero(vertexMass.size());
//
//	for (size_t iE = 0; iE < nEdges; iE++)
//	{
//		FloatingType edgeLength = orgLengths(iE);
//		vertexAdjacentELengthsAll(edges(0, iE)) += edgeLength / 2;
//		vertexAdjacentELengthsAll(edges(1, iE)) += edgeLength / 2;
//	}
//
//	for (size_t iV = 0; iV < numVertices(); iV++)
//	{
//		vertexMass(iV) = totalMass * vertexAdjacentELengthsAll(iV) / totalSpringLen;
//	}
//
//	assert(abs(totalMass - (vertexMass).sum()) < 1e-4);
//	vertexInvMass.array() = 1.f / vertexMass.array();
//
//	for (size_t iV = 0; iV < numVertices(); iV++)
//	{
//		vertextInvMassX3(iV * 3) = vertexInvMass(iV);
//		vertextInvMassX3(iV * 3 + 1) = vertexInvMass(iV);
//		vertextInvMassX3(iV * 3 + 2) = vertexInvMass(iV);
//
//		vertexMassX3(iV * 3) = vertexMass(iV);
//		vertexMassX3(iV * 3 + 1) = vertexMass(iV);
//		vertexMassX3(iV * 3 + 2) = vertexMass(iV);
//	}
//}
//
//
//void GAIA::EBDTetMesh_MassSpring::forwardStepForwardEuler()
//{
//	positionsPrevFlatten = positionsFlatten;
//	velocityPrevFlatten = velocityFlatten;
//
//	// std::cout << "internalForceFlatten" << internalForceFlatten.transpose() << "\n";
//	positionsFlatten += velocityFlatten * pPhysicsParams->dt;
//
//	velocityFlatten += pPhysicsParams->dt * internalForceFlatten.cwiseProduct(vertextInvMassX3);
//	handleVelConstraints();
//
//	
//}
//
//void GAIA::EBDTetMesh_MassSpring::evaluateInternalForce(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> internalForceOut, bool addGravity)
//{
//	internalForceOut.setZero();
//
//	for (size_t iEdge = 0; iEdge < nEdges; iEdge++)
//	{
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = position.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - position.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//
//		FloatingType l0 = orgLengths(iEdge);
//
//		Vec3 force = (adjustedSpringStiffness(iEdge) * (l0 - l) / l) * diff;
//
//		internalForceOut.block<3, 1>(v1 * POINT_VEC_DIMS, 0) += force;
//		internalForceOut.block<3, 1>(v2 * POINT_VEC_DIMS, 0) -= force;
//
//	}
//	if (addGravity) {
//		for (size_t iVert = 0; iVert < numVertices(); iVert++)
//		{
//			internalForceOut.block<3, 1>(iVert * POINT_VEC_DIMS, 0) += pPhysicsParams->gravity * vertexMass(iVert);
//		}
//	}
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateInternalPotentialEnergy(Eigen::Ref<VecDynamic> position)
//{
//	return evaluateSpringEnergy(position);
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateInternalPotentialEnergy(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> perElementEnergy)
//{
//	FloatingType e = 0.f;
//	for (size_t iEdge = 0; iEdge < nEdges; iEdge++)
//	{
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = position.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - position.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		perElementEnergy(iEdge) = 0.5f * adjustedSpringStiffness(iEdge) * (l - l0) * (l - l0);
//
//		e += perElementEnergy(iEdge);
//	}
//	return e;
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateMeritEnergyPervertex(int iV)
//{
//	FloatingType inertiaEnergy = vertexMass(iV) * (vertex(iV) - BLOCK_DIM3x1_AT(inertia, iV)).squaredNorm() / (2 * pPhysicsParams->dt * pPhysicsParams->dt);
//	return inertiaEnergy + evaluateLocalInternalEnergyPerVertex(iV);
//}
//
//void GAIA::EBDTetMesh_MassSpring::evaluateSpringEnergyPerSpring(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> springEnergy)
//{
//	cpu_parallel_for(0, nEdges, [&](int iEdge) {
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = position.block<3,1>(v1*POINT_VEC_DIMS, 0) - position.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		springEnergy(iEdge)  = 0.5f * adjustedSpringStiffness(iEdge) * (l - l0) * (l - l0);
//	});
//
//
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateSpringEnergy(Eigen::Ref<VecDynamic> position)
//{
//	FloatingType e = 0.f;
//	for (size_t iEdge = 0; iEdge < nEdges; iEdge++)
//	{
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = position.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - position.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		e += 0.5 * adjustedSpringStiffness(iEdge) * (l - l0) * (l - l0);
//	}
//	return e;
//}
//
//
//void GAIA::EBDTetMesh_MassSpring::computeElasiticForce(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4, const Mat3& DsInv, 
//	FloatingType restVol, Vec3& p1Grad, Vec3& p2Grad, Vec3& p3Grad, Vec3& p4Grad)
//{
//
//}
//
//
//
//GAIA::FloatingType GAIA::EBDTetMesh_MassSpring::evaluateEnergy(int iTet)
//{
//	return 0;
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateElasticEnergy()
//{
//
//	return evaluateInternalPotentialEnergy(positionsFlatten);
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::adjustedSpringStiffness(int iSpring)
//{
//	//return pObjectParamsMaterial->springStiffness / (orgLengths(iSpring) / avgSpringLength);
//	return orgLengths(iSpring) * pObjectParamsMaterial->springStiffness / avgSpringLength;
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateLocalInternalEnergyPerVertex(int vId)
//{
//	FloatingType energy = 0.f;
//	const std::vector<int>& neiEdges = vertexNeighborEdges[vId];
//	for (size_t iNeiE = 0; iNeiE < neiEdges.size(); iNeiE++)
//	{
//		int iEdge = neiEdges[iNeiE];
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = positionsFlatten.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - positionsFlatten.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		energy += 0.5 * adjustedSpringStiffness(iEdge) * (l - l0) * (l - l0);
//	}
//	return energy;
//}
//
//void GAIA::EBDTetMesh_MassSpring::evaluateLocalInternalEnergyHessianPerVertex(int vId, Mat3 & h)
//{
//	FloatingType energy = 0.f;
//	const std::vector<int>& neiEdges = vertexNeighborEdges[vId];
//	h = Mat3::Zero();
//	for (size_t iNeiE = 0; iNeiE < neiEdges.size(); iNeiE++)
//	{
//		int iEdge = neiEdges[iNeiE];
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = positionsFlatten.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - positionsFlatten.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//		// evaluate hessian
//		Mat3 h_1_1;
//		h_1_1 = adjustedSpringStiffness(iEdge) * (Mat3::Identity() - (l0 / l) * (Mat3::Identity() - diff * diff.transpose() / (l * l)));
//		h += h_1_1;
//
//	}
//	return ;
//}
//
//FloatingType GAIA::EBDTetMesh_MassSpring::evaluateLocalInternalEnergyPerVertex(int vId, const Vec3& newPos, const Eigen::Ref<VecDynamic> allVerts)
//{
//	FloatingType energy = 0.f;
//	const std::vector<int>& neiVerts = vertexNeighborEdges[vId];
//	for (size_t iNeiE = 0; iNeiE < neiVerts.size(); iNeiE++)
//	{
//		int iEdge = neiVerts[iNeiE];
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff;
//		if (v1 == vId)
//		{
//			diff = newPos - allVerts.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		}
//		else {
//			diff = allVerts.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - newPos;
//		}
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		energy += 0.5 * adjustedSpringStiffness(iEdge) * (l - l0) * (l - l0);
//	}
//	return energy;
//}
//
//void GAIA::EBDTetMesh_MassSpring::evaluateInternalForcePerVertex(int vId, Vec3& force, bool addGravity)
//{
//	force.setZero();
//	const std::vector<int>& neiVerts = vertexNeighborEdges[vId];
//	for (size_t iNeiE = 0; iNeiE < neiVerts.size(); iNeiE++)
//	{
//		int iEdge = neiVerts[iNeiE];
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//		Vec3 diff = positionsFlatten.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - positionsFlatten.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//		FloatingType l = diff.norm();
//		FloatingType l0 = orgLengths(iEdge);
//
//		if (v1 == vId)
//		{
//			// force is the negative of energy gradient
//			force += (adjustedSpringStiffness(iEdge) * (l0 - l) / l) * diff;
//
//		}
//		else {
//			// force is the negative of energy gradient 
//			force -= (adjustedSpringStiffness(iEdge) * (l0 - l) / l) * diff;
//		}
//	}
//	if (addGravity)
//	{
//		force += pPhysicsParams->gravity * vertexMass(vId);
//	}
//}
//
//void GAIA::EBDTetMesh_MassSpring::evaluateInternalForcePerVertex(int vId, const Vec3& newPos, const Eigen::Ref<VecDynamic> allVerts, Vec3& force, bool addGravity)
//{
//	force.setZero();
//	const std::vector<int>& neiVerts = vertexNeighborEdges[vId];
//	for (size_t iNeiE = 0; iNeiE < neiVerts.size(); iNeiE++)
//	{
//		int iEdge = neiVerts[iNeiE];
//		IdType v1 = edges(0, iEdge);
//		IdType v2 = edges(1, iEdge);
//
//		if (v1 == vId)
//		{
//			Vec3 diff = newPos - allVerts.block<3, 1>(v2 * POINT_VEC_DIMS, 0);
//			FloatingType l = diff.norm();
//			FloatingType l0 = orgLengths(iEdge);
//			force += (adjustedSpringStiffness(iEdge) * (l0 - l) / l) * diff;
//		}
//		else {
//			Vec3 diff = allVerts.block<3, 1>(v1 * POINT_VEC_DIMS, 0) - newPos;
//			FloatingType l = diff.norm();
//			FloatingType l0 = orgLengths(iEdge);
//			force += (adjustedSpringStiffness(iEdge) * (l0 - l) / l) * diff;
//		}
//	}
//	if (addGravity)
//	{
//	force += pPhysicsParams->gravity * vertexMass(vId);
//	}
//
//}
//
//bool GAIA::ObjectParametersEBDMassSpring::fromJson(nlohmann::json& objectJsonParams)
//{
//	ObjectParamsVBD::fromJson(objectJsonParams);
//	EXTRACT_FROM_JSON(objectJsonParams, springStiffness);
//	EXTRACT_FROM_JSON(objectJsonParams, damping);
//	EXTRACT_FROM_JSON(objectJsonParams, lengthBasedMass);
//	EXTRACT_FROM_JSON(objectJsonParams, localDescent);
//
//	return true;
//}
//
//bool GAIA::ObjectParametersEBDMassSpring::toJson(nlohmann::json& objectJsonParams)
//{
//	ObjectParamsVBD::toJson(objectJsonParams);
//	PUT_TO_JSON(objectJsonParams, springStiffness);
//	PUT_TO_JSON(objectJsonParams, damping);
//	PUT_TO_JSON(objectJsonParams, lengthBasedMass);
//	PUT_TO_JSON(objectJsonParams, localDescent);
//
//	return true;
//}
