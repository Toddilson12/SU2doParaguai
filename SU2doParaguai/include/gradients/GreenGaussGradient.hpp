#pragma once

#include "Gradient.hpp"

namespace FVGradient
{
	class CGreenGaussGradient : public CGradient //DONT USE IN MIXED GRIDS!!!
	{
	public:
		CGreenGaussGradient(FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
			: CGradient(mesh, variables) { }
		void ComputeGradients(const GradOptions& opts, ulong iVarBegin = 0, ulong iVarEnd = -1) override;
	};


}

#include "GreenGaussGradient.inl"