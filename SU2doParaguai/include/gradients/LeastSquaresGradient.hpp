#pragma once
#include "Gradient.hpp"

namespace FVGradient
{
	class CLeastSquaresGradient : public CGradient
	{
	private:
		// Inverse matrix
		FVUtil::Block2DMatrix invKMat; // this is the A^T W^2 A. Only stores 3 elements. Order is: (1,1), (1,2), (2,2) 
		//I could technically store the weights (not W) for each deltaPhi, but that would be a 2xm matrix for each node
		// where m is the num of neighbors. This costs more memory (bottleneck)
		//std::vector<ulong> illNodes; // ill-conditioned nodes, marked for QR decomp
		double* weights; //equivalent to w_ij^2 = 1 / (delta_x^2 + delta_y^2)
	public:
		CLeastSquaresGradient(const GradOptions& opts, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables);
		void ComputeGradients(const GradOptions& opts, ulong iVarBegin = -1, ulong iVarEnd = -1);
		~CLeastSquaresGradient()
		{
			invKMat.DeallocateMemory();
			delete[] weights;
		}
	};

}
#include "LeastSquaresGradient.inl"