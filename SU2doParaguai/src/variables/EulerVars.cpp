#include "../../include/variables/EulerVars.hpp"

namespace FVVariables
{
	CEulerVariable::CEulerVariable(FVSolver::CConfig* config, ulong n_poin)
		: CVariableContainer(config, n_poin, 2, 4, 4, 4)
	{
		//Note: CIndices is already init to correct idxs by default for Euler variables
		
		//Initialize all matrices
		for (ulong i = 0; i < nPoin; i++)
		{
			for (ulong j = 0; j < nVars; j++) // note: nVars = nPrims = nPrimsGrad
			{
				for (ulong d = 0; d < nDim; d++)
				{
					SetGradient(i, j, static_cast<FVUtil::COORDINATES>(d), 0.0);
					SetGradientPrim(i, j, static_cast<FVUtil::COORDINATES>(d), 0.0);
				}
				SetSolution(i, j, 0.0);
				SetSolutionOld(i, j, 0.0);
				SetPrimitive(i, j, 0.0);
				SetResidual(i, j, 0.0);
				SetSource(i, j, 0.0);
				SetSlopeLimiter(i, j, 0.0);
			}
		}
		
	}

	CEulerVariable::~CEulerVariable()
	{
		//delete any new pointers
	}


}