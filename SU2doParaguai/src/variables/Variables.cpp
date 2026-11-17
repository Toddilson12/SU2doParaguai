#include "../../include/variables/Variables.hpp"


namespace FVVariables
{
	CVariableContainer::CVariableContainer(FVSolver::CConfig* config, ulong n_poin, int n_dim, ulong n_vars, ulong n_prims, ulong n_primsgrad)
		: nPoin(n_poin), nDim(n_dim), nVars(n_vars), nPrims(n_prims), nPrimsGrad(n_primsgrad), Indices(CIndices()) //necessary? no...
	{
		//allocate memory for the fields
		SolutionMat = Matrix::AllocateBlock2DMatrix(nPoin, nVars);
		SolutionOldMat = Matrix::AllocateBlock2DMatrix(nPoin, nVars);

		PrimitiveMat = Matrix::AllocateBlock2DMatrix(nPoin, nPrims);
		ResidualMat = Matrix::AllocateBlock2DMatrix(nPoin, nVars);

		GradientMat = Matrix::AllocateBlock2DMatrix(nPoin, nVars * nDim);
		GradientPrimMat = Matrix::AllocateBlock2DMatrix(nPoin, nPrimsGrad * nDim);
	
		SourceMat = Matrix::AllocateBlock2DMatrix(nPoin, nVars);

		SlopeLimiterMat = Matrix::AllocateBlock2DMatrix(nPoin, nPrims);
	}

	CVariableContainer::~CVariableContainer()
	{
		SolutionMat.DeallocateMemory();
		SolutionOldMat.DeallocateMemory();
		PrimitiveMat.DeallocateMemory();
		ResidualMat.DeallocateMemory();
		GradientMat.DeallocateMemory();
		GradientPrimMat.DeallocateMemory();
		SourceMat.DeallocateMemory();
		SlopeLimiterMat.DeallocateMemory();
	}
	

}
