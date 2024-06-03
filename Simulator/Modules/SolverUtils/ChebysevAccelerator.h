#pragma once
#include "../Types/Types.h"

namespace GAIA
{
	struct ChebysevAccelerator {
		typedef std::shared_ptr<ChebysevAccelerator> SharedPtr;
		typedef ChebysevAccelerator* Ptr;

		ChebysevAccelerator(FloatingType rho_in)
			: rho(rho_in)
		{}

		// order starts from 1
		FloatingType getAcceleratorOmega(int order, CFloatingType prevOmega)
		{
			if (rho==0.f)
			{
				return 1.f;
			}
			switch (order)
			{
			case 1:
				return 1.f;
				break;
			case 2:
				return  2.f / (2.f - SQR(rho));
				break;
			default:
				assert(order > 0);

				return 4.f / (4.f - SQR(rho) * prevOmega);;
				break;
			}
		}

		CFloatingType rho;
	};
}
