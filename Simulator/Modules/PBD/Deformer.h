#pragma once 

namespace GAIA {
	class PBDPhysics;

	struct BaseDeformer
	{
	public:
		BaseDeformer() {};
		~BaseDeformer() {};

		virtual void operator()(PBDPhysics* physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) = 0;

	private:

	};

	typedef std::shared_ptr<BaseDeformer> DeformerPtr;
}