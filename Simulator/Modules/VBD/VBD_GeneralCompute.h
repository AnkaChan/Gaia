#pragma once
#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrix.h"
#include "VBD_BaseMeshGPU.h"

namespace GAIA {

	GPU_CPU_INLINE_FUNC  void accumulateInertiaForceAndHessian(VBDBaseTetMeshGPU* pTetMeshGPU, int iV, CFloatingTypeGPU dtSqrReciprocal,
		FloatingTypeGPU* force, FloatingTypeGPU* h);

    GPU_CPU_INLINE_FUNC  void accumulateVertexFrictionGPU(CFloatingTypeGPU mu, CFloatingTypeGPU lambda, CFloatingTypeGPU* T_3x2, CFloatingTypeGPU* u_2x1,
        CFloatingTypeGPU epsU, FloatingTypeGPU* force, FloatingTypeGPU* hessian)
    {
        // Friction
        CFloatingTypeGPU uNorm = CuMatrix::vec2Norm(u_2x1);
        if (uNorm > 0)
        {
            // IPC friction 
            // https://github.com/ipc-sim/ipc-toolkit/blob/main/src/ipc/friction/smooth_friction_mollifier.cpp
            FloatingTypeGPU f1_SF_over_x;
            FloatingTypeGPU df1_x_minus_f1_over_x3;
            if (uNorm > epsU)
            {
                f1_SF_over_x = 1 / uNorm;
                df1_x_minus_f1_over_x3 = -1 / (uNorm * uNorm * uNorm);
            }
            else
            {
                f1_SF_over_x = (-uNorm / epsU + 2) / epsU;
                df1_x_minus_f1_over_x3 = -1 / (uNorm * epsU * epsU);
            }
            FloatingTypeGPU forceFriction[3];
            CuMatrix::matMulMxN<FloatingTypeGPU, 3, 2, 1>(T_3x2, u_2x1, forceFriction);

            // force = -mu * lambda * f1_SF_over_x * T_3x2 * u_2x1;
            CuMatrix::vec3MulAddTo(forceFriction, -mu * lambda * f1_SF_over_x, force);

            FloatingTypeGPU f1_SF_over_x_I2[4] = {
                mu * lambda * f1_SF_over_x, // (0,0)
                0.f, // (1, 0)
                0.f, // (0, 1)
                mu * lambda * f1_SF_over_x // (1, 1)
            };

            // compute T * (f1_SF_over_x * Mat2::Identity()) * T.transpose();
            FloatingTypeGPU m1[6];
            CuMatrix::matMulMxN<FloatingTypeGPU, 3, 2, 2>(T_3x2, f1_SF_over_x_I2, m1);

            for (int r = 0; r < 3; r++)
            {
                for (int c = 0; c < 3; c++) {

                    for (int i = 0; i < 2; i++) {
                        CuMatrix::accessMatElement<FloatingTypeGPU, 3, 3>(hessian, r, c) +=
                            CuMatrix::accessMatElement<FloatingTypeGPU, 3, 2>(m1, r, i)
                            * CuMatrix::accessMatElement<FloatingTypeGPU, 3, 2>(T_3x2, c, i);
                        //  CuMatrix::accessMatElement<FloatingTypeGPU, 3, 2>(T_3x2, c, i) 
                        //  == CuMatrix::accessMatElement<FloatingTypeGPU, 2, 3>(T_3x2_t, i, c) 
                    }
                }
            }

            //hessian = mu * lambda * T * (df1_x_minus_f1_over_x3 * u * u.transpose() + f1_SF_over_x * Mat2::Identity()) * T.transpose();
        }

        return;
    }

	GPU_CPU_INLINE_FUNC  void GAIA::accumulateInertiaForceAndHessian(VBDBaseTetMeshGPU* pTetMeshGPU, int iV, CFloatingTypeGPU dtSqrReciprocal,
		FloatingTypeGPU* force, FloatingTypeGPU* h)
	{
		CFloatingTypeGPU vertexMass = pTetMeshGPU->vertexMass[iV];
		CFloatingTypeGPU* inertia = pTetMeshGPU->inertia + 3 * iV;
		CFloatingTypeGPU* pos = pTetMeshGPU->vertPos + 3 * iV;

		FloatingTypeGPU inertiaForce[3];
		CFloatingTypeGPU vertexMass_X_dtSqrReciprocal = dtSqrReciprocal * vertexMass;

		CuMatrix::vec3Minus(inertia, pos, inertiaForce);
		CuMatrix::vec3Mul(inertiaForce, vertexMass_X_dtSqrReciprocal, inertiaForce);
		CuMatrix::vec3Add(force, inertiaForce, force);

		h[0] += vertexMass_X_dtSqrReciprocal;
		h[4] += vertexMass_X_dtSqrReciprocal;
		h[8] += vertexMass_X_dtSqrReciprocal;
	}

	GPU_CPU_INLINE_FUNC void accumulateBoundaryForceAndHessianGPU(VBDPhysicsDataGPU* pPhysicsData, VBDBaseTetMeshGPU* pTetMeshGPU, int iV,
		FloatingTypeGPU* force,	FloatingTypeGPU* hessian) {
		CFloatingTypeGPU boundaryCollisionStiffness = pPhysicsData->boundaryCollisionStiffness;
        bool frictionApplied = false;
		for (size_t iDim = 0; iDim < 3; iDim++)
		{
			// FloatingTypeGPU  contactNormal[3] = { 0.f, 0.f, 0.f };

			CFloatingTypeGPU lowerBound = pPhysicsData->worldBounds[iDim];
			CFloatingTypeGPU upperBound = pPhysicsData->worldBounds[iDim+3];

			CFloatingTypeGPU* v = pTetMeshGPU->vertPos + VERTEX_BUFFER_STRIDE * iV;
			if (v[iDim] < lowerBound)
			{
                // Penalty Force
                CFloatingTypeGPU penetrationDepth = lowerBound - v[iDim];
				force[iDim] += penetrationDepth * boundaryCollisionStiffness;
				hessian[iDim + iDim * 3] += boundaryCollisionStiffness;
                // Friction
                if (!frictionApplied)
                    // avoid multiple friction when there it exeeds more than one boundary
                {
                    frictionApplied = true;
                    FloatingTypeGPU dx[3]; 
                    CuMatrix::vec3Minus(pTetMeshGPU->getVert(iV), pTetMeshGPU->getVertPrevPos(iV), dx);

                    FloatingTypeGPU T[6] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

                    T[(iDim + 1) % 3] = 1.f;
                    T[(iDim + 2) % 3 + 3] = 1.f;

                    FloatingTypeGPU u[2] = {
                        T[0] * dx[0] + T[1] * dx[1] + T[2] * dx[2],
                        T[3] * dx[0] + T[4] * dx[1] + T[5] * dx[2],
                    };

                    CFloatingTypeGPU lambda = penetrationDepth * boundaryCollisionStiffness;
                    CFloatingTypeGPU mu = pPhysicsData->boundaryFrictionDynamic;
                    CFloatingTypeGPU epsV = pPhysicsData->boundaryFrictionEpsV;
                    CFloatingTypeGPU epsU = epsV * pPhysicsData->dt;

                    accumulateVertexFrictionGPU(mu, lambda, T, u, epsU, force, hessian);
                }

			}
			else if (v[iDim] > upperBound)
			{
				CFloatingTypeGPU penetrationDepth = v[iDim] - upperBound;
				force[iDim] -= penetrationDepth * boundaryCollisionStiffness;

				hessian[iDim + iDim * 3] += boundaryCollisionStiffness;

                // Friction
                if (!frictionApplied)
                // avoid multiple friction when there it exeeds more than one boundary
                {
                    frictionApplied = true;
                    FloatingTypeGPU dx[3];
                    CuMatrix::vec3Minus(pTetMeshGPU->getVert(iV), pTetMeshGPU->getVertPrevPos(iV), dx);

                    FloatingTypeGPU T[6] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

                    T[(iDim + 1) % 3] = 1.f;
                    T[(iDim + 2) % 3 + 3] = 1.f;

                    FloatingTypeGPU u[2] = {
                        T[0] * dx[0] + T[1] * dx[1] + T[2] * dx[2],
                        T[3] * dx[0] + T[4] * dx[1] + T[5] * dx[2],
                    };

                    CFloatingTypeGPU lambda = penetrationDepth * boundaryCollisionStiffness;
                    CFloatingTypeGPU mu = pPhysicsData->boundaryFrictionDynamic;
                    CFloatingTypeGPU epsV = pPhysicsData->boundaryFrictionEpsV;
                    CFloatingTypeGPU epsU = epsV * pPhysicsData->dt;

                    accumulateVertexFrictionGPU(mu, lambda, T, u, epsU, force, hessian);
                }
            }
		}
	}

    enum class ClosestPointOnTriangleTypeGPU
    {
        AtA,
        AtB,
        AtC,
        AtAB,
        AtBC,
        AtAC,
        AtInterior,
        NotFound
    };

    GPU_CPU_INLINE_FUNC ClosestPointOnTriangleTypeGPU closestPointTriangle(CFloatingTypeGPU * p,
        CFloatingTypeGPU* a, CFloatingTypeGPU* b, CFloatingTypeGPU* c, 
        FloatingTypeGPU* baryCentrics, FloatingTypeGPU* closestPt)
    {
        ClosestPointOnTriangleTypeGPU pointType;
        FloatingTypeGPU ab[3]; //b - a;
        CuMatrix::vec3Minus(b, a, ab);
        FloatingTypeGPU ac[3]; //c - a;
        CuMatrix::vec3Minus(c, a, ac);
        FloatingTypeGPU ap[3]; //p - a;
        CuMatrix::vec3Minus(p, a, ap);

        CFloatingTypeGPU d1 = CuMatrix::vec3DotProduct(ab, ap);
        CFloatingTypeGPU d2 = CuMatrix::vec3DotProduct(ac, ap);
        if (d1 <= 0.f && d2 <= 0.f) {
            pointType = ClosestPointOnTriangleTypeGPU::AtA;
            CuMatrix::vec3Set(baryCentrics, 1.f, 0.f, 0.f);
            CuMatrix::vec3Set(closestPt, a);
            return pointType;
        }

        FloatingTypeGPU bp[3]; // = p - b;
        CuMatrix::vec3Minus(p, b, bp);
        CFloatingTypeGPU d3 = CuMatrix::vec3DotProduct(ab, bp);
        CFloatingTypeGPU d4 = CuMatrix::vec3DotProduct(ac, bp);
        if (d3 >= 0.f && d4 <= d3) {
            pointType = ClosestPointOnTriangleTypeGPU::AtB;
            CuMatrix::vec3Set(baryCentrics, 0.f, 1.f, 0.f);
            CuMatrix::vec3Set(closestPt, b);
            return pointType;
        }

        FloatingTypeGPU cp[3]; // = p - c;
        CuMatrix::vec3Minus(p, c, cp);
        CFloatingTypeGPU d5 = CuMatrix::vec3DotProduct(ab, cp);
        CFloatingTypeGPU d6 = CuMatrix::vec3DotProduct(ac, cp);
        if (d6 >= 0.f && d5 <= d6) {
            pointType = ClosestPointOnTriangleTypeGPU::AtC;
            CuMatrix::vec3Set(baryCentrics, 0.f, 0.f, 1.f);
            CuMatrix::vec3Set(closestPt, c);
            return pointType;
        }

        CFloatingTypeGPU vc = d1 * d4 - d3 * d2;
        if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
        {
            CFloatingTypeGPU v = d1 / (d1 - d3);
            pointType = ClosestPointOnTriangleTypeGPU::AtAB;
            CuMatrix::vec3Set(baryCentrics, 1.0f - v, v, 0.f);
            CuMatrix::vec3lerp(a, v, ab, closestPt);
            //return a + v * ab;
            return pointType;
        }

        CFloatingTypeGPU vb = d5 * d2 - d1 * d6;
        if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
        {
            CFloatingTypeGPU v = d2 / (d2 - d6);
            pointType = ClosestPointOnTriangleTypeGPU::AtAC;
            CuMatrix::vec3Set(baryCentrics, 1.0f - v, 0.f, v);
            CuMatrix::vec3lerp(a, v, ac, closestPt);
            //return a + v * ac;
            return pointType;
        }

        CFloatingTypeGPU va = d3 * d6 - d5 * d4;
        if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
        {
            pointType = ClosestPointOnTriangleTypeGPU::AtBC;
            CFloatingTypeGPU v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            CuMatrix::vec3Set(baryCentrics, 0.f, 1.f - v, v);
            //return b + v * (c - b);
            FloatingTypeGPU bc[3]; // = p - c;
            CuMatrix::vec3Minus(c, b, bc);
            CuMatrix::vec3lerp(b, v, bc, closestPt);
            return pointType;
        }

        CFloatingTypeGPU denom = 1.f / (va + vb + vc);
        CFloatingTypeGPU v = vb * denom;
        CFloatingTypeGPU w = vc * denom;
        pointType = ClosestPointOnTriangleTypeGPU::AtInterior;
        CuMatrix::vec3Set(baryCentrics, 1.f - v - w, v, w);
        CuMatrix::vec3lerp(a, v, ab, closestPt);
        CuMatrix::vec3MulAddTo(ac, w, closestPt);

        // return a + v * ab + w * ac;
        return pointType;
    }


}