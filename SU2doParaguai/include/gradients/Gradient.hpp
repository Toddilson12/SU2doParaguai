#pragma once

#include "../common.hpp"
#include "../variables/Variables.hpp"
#include "../mesh/mesh.hpp"

namespace FVGradient
{
	struct GradOptions
	{
		bool weighted = false; // if true, then w_ij = 1/d_ij => w_ij^2 = 1/... with no sqrt!!! If false, then the whole thing is just one
		bool computePrim = true;
	};

	class CGradient
	{
	protected:
		FVVariables::CVariableContainer* Variables;
		FVMesh::CMesh* Mesh;
	public:
		CGradient(FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables);
		virtual ~CGradient();
		
		virtual void ComputeGradients(const GradOptions& opts, ulong iVarBegin = -1, ulong iVarEnd = -1) = 0;
	};


}