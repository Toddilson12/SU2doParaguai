#include "../../include/gradients/Gradient.hpp"


namespace FVGradient
{
	CGradient::CGradient(FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
	{
		this->Mesh = mesh;
		this->Variables = variables;
	}
	CGradient::~CGradient()
	{
		this->Mesh = nullptr;
		this->Variables = nullptr;
	}

}
