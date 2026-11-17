#pragma once
#include "../common.hpp"
#include "../config/config.hpp"
#include "Indices.hpp"


namespace FVVariables
{
	class CVariableContainer
	{
	protected:
		using Matrix = FVUtil::Block2DMatrix;
		//solution (usually conservative vars)
		Matrix SolutionMat;
		Matrix SolutionOldMat; //time level n-1
		//primitives
		Matrix PrimitiveMat;
		//residual
		Matrix ResidualMat;
		//gradient for solution
		Matrix GradientMat;
		//gradient for prim
		Matrix GradientPrimMat;
		//source terms
		Matrix SourceMat;
		//limiters
		Matrix SlopeLimiterMat;
		//indices (its in the base class so that it can be referred to in other classes)
		//maybe add something for cfl/time stepping? maybe conv wave speed?

		//sizes
		ulong nPoin = 0;
		int nDim = 0;

		ulong nVars = 0;
		ulong nPrims = 0;

		ulong nPrimsGrad = 0; // only compute and store some of the prim grads (for slope limiters)

		
	public:
		CIndices Indices;
		CVariableContainer(FVSolver::CConfig* config, ulong n_poin, int n_dim, ulong n_vars, ulong n_prims, ulong n_primsgrad);
		virtual ~CVariableContainer();
		//dont forget to inline
		// 
		// Setters (global and at node)

		ulong GetNumVars()
		{
			return nVars;
		}

		ulong GetNumPrims()
		{
			return nPrims;
		}

		ulong GetNumPrimsGrad()
		{
			return nPrimsGrad;
		}



		void SetSolution(double* solution) // just copying... no ownership change
		{
			SolutionMat = solution;
		}
		void SetSolution(ulong i, double* solution_vec)
		{
			for (ulong j = 0; j < nVars; j++)
				SolutionMat(i, j) = solution_vec[j];
		}
		void SetSolution(ulong i, ulong iVar, double solution_val)
		{
			SolutionMat(i, iVar) = solution_val;
		}

		void SetSolutionOld(double* solution) // just copying... no ownership change
		{
			SolutionOldMat = solution;
		}
		void SetSolutionOld(ulong i, double* solution_vec)
		{
			for (ulong j = 0; j < nVars; j++)
				SolutionOldMat(i, j) = solution_vec[j];
		}
		void SetSolutionOld(ulong i, ulong iVar, double solution_val)
		{
			SolutionOldMat(i, iVar) = solution_val;
		}


		void SetPrimitive(double* primitives)
		{
			PrimitiveMat = primitives;
		}
		void SetPrimitive(ulong i, double* primitive_vec)
		{
			for (ulong j = 0; j < nPrims; j++)
				PrimitiveMat(i, j) = primitive_vec[j];
		}
		void SetPrimitive(ulong i, ulong iVar, double primitive_val)
		{
			PrimitiveMat(i, iVar) = primitive_val;
		}


		void SetResidual(double* residuals)
		{
			ResidualMat = residuals;
		}
		void SetResidual(ulong i, double* residual_vec)
		{
			for (ulong j = 0; j < nVars; j++)
				ResidualMat(i, j) = residual_vec[j];
		}
		void SetResidual(ulong i, ulong iVar, double residual_val)
		{
			ResidualMat(i, iVar) = residual_val;
		}

		void SetGradient(double* gradients)
		{
			GradientMat = gradients;
		}


		void SetGradient(ulong i, double* gradient_vec) //DO NOT USE
		{
			for (ulong j = 0; j < nVars; j++)
				GradientMat(i, j) = gradient_vec[j];
		}

		void SetGradient(ulong i, ulong raw_idx, double gradient_val) // DO NOT USE
		{
			GradientMat(i, raw_idx) = gradient_val;
		}

		void SetGradient(ulong i, ulong iVar, FVUtil::COORDINATES dir, double gradient_val)
		{
			GradientMat(i, iVar * (ulong)nDim + (ulong)dir) = gradient_val;
		}
		
		void SetGradient(ulong i, ulong iVar, double* gradient_vec) //No length check
		{
			for (ulong iDim = 0; iDim < nDim; iDim++)
				SetGradient(i, iVar, static_cast<FVUtil::COORDINATES>(iDim), gradient_vec[iDim]); //optmize later
		}

		void SetGradient(ulong i, ulong iVar, const FVUtil::Vec2& gradient_vec) //ignore third component
		{
			SetGradient(i, iVar, FVUtil::COORDINATES::X, gradient_vec.x);
			SetGradient(i, iVar, FVUtil::COORDINATES::Y, gradient_vec.y);
		}

		void AddGradient(ulong i, ulong iVar, const FVUtil::Vec2& gradient_vec)
		{
			GradientMat(i, iVar * (ulong)nDim + (ulong)FVUtil::COORDINATES::X) += gradient_vec.x;
			GradientMat(i, iVar * (ulong)nDim + (ulong)FVUtil::COORDINATES::Y) += gradient_vec.y;
		}

		void ScaleGradient(ulong i, ulong iVar, double scaling_factor)
		{
			for(ulong iDim = 0; iDim < nDim; iDim++)
				GradientMat(i, iVar * (ulong)nDim + iDim) *= scaling_factor;
		}
		//talvez uma das coisas mais idiotas que eu ja escrevi

		void SetGradientPrim(double* gradients)
		{
			GradientPrimMat = gradients;
		}

		void SetGradientPrim(ulong i, double* gradient_vec)
		{
			for (ulong j = 0; j < nPrimsGrad; j++)
				GradientPrimMat(i, j) = gradient_vec[j];
		}

		void SetGradientPrim(ulong i, ulong raw_idx, double gradient_val)
		{
			GradientPrimMat(i, raw_idx) = gradient_val;
		}

		void SetGradientPrim(ulong i, ulong iPrim, FVUtil::COORDINATES dir, double gradient_val)
		{
			GradientPrimMat(i, iPrim * (ulong)nDim + (ulong)dir) = gradient_val;
		}

		void SetGradientPrim(ulong i, ulong iPrim, double* gradient_vec) //No length check
		{
			for (ulong iDim = 0; iDim < nDim; iDim++)
				SetGradientPrim(i, iPrim, static_cast<FVUtil::COORDINATES>(iDim), gradient_vec[iDim]); //optmize later
		}

		void SetGradientPrim(ulong i, ulong iPrim, const FVUtil::Vec2& gradient_vec) //ignore third component
		{
			SetGradientPrim(i, iPrim, FVUtil::COORDINATES::X, gradient_vec.x);
			SetGradientPrim(i, iPrim, FVUtil::COORDINATES::Y, gradient_vec.y);
		}

		void AddGradientPrim(ulong i, ulong iPrim, const FVUtil::Vec2& gradient_vec)
		{
			GradientPrimMat(i, iPrim * (ulong)nDim + (ulong)FVUtil::COORDINATES::X) += gradient_vec.x;
			GradientPrimMat(i, iPrim * (ulong)nDim + (ulong)FVUtil::COORDINATES::Y) += gradient_vec.y;
		}

		void ScaleGradientPrim(ulong i, ulong iPrim, double scaling_factor)
		{
			for(ulong iDim = 0; iDim < nDim; iDim++)
				GradientPrimMat(i, iPrim * (ulong)nDim + iDim) *= scaling_factor;
		}


		void SetSource(double* sources)
		{
			SourceMat = sources;
		}

		void SetSource(ulong i, double* source_vec)
		{
			for (ulong j = 0; j < nVars; j++)
				SourceMat(i, j) = source_vec[j];
		}

		void SetSource(ulong i, ulong iVar, double source_val)
		{
			SourceMat(i, iVar) = source_val;
		}

		void SetSlopeLimiter(double* slope_limiter)
		{
			SlopeLimiterMat = slope_limiter;
		}

		void SetSlopeLimiter(ulong i, double* slope_vec)
		{
			for (ulong j = 0; j < nPrims; j++)
				SlopeLimiterMat(i, j) = slope_vec[j];
		}

		void SetSlopeLimiter(ulong i, ulong iVar, double slope_val)
		{
			SlopeLimiterMat(i, iVar) = slope_val;
		}


		// Getters (global and at node)
		
		double* GetSolution()
		{
			return SolutionMat.block_ptr;
		}

		double* GetSolution(ulong i)
		{
			return SolutionMat(i);
		}

		double GetSolution(ulong i, ulong iVar)
		{
			return SolutionMat(i, iVar);
		}

		double* GetPrimitive()
		{
			return PrimitiveMat.block_ptr;
		}
		double* GetPrimitive(ulong i)
		{
			return PrimitiveMat(i);
		}
		double GetPrimitive(ulong i, ulong iVar)
		{
			return PrimitiveMat(i, iVar);
		}
		
		double* GetResidual()
		{
			return ResidualMat.block_ptr;
		}
		double* GetResidual(ulong i)
		{
			return ResidualMat(i);
		}
		double GetResidual(ulong i, ulong iVar)
		{
			return ResidualMat(i, iVar);
		}

		double* GetGradient()
		{
			return GradientMat.block_ptr;
		}
		double* GetGradient(ulong i)
		{
			return GradientMat(i);
		}
		double GetGradient(ulong i, ulong iVar) // NAO FUNCIONA MAIS!!!
		{
			return GradientMat(i, iVar);
		}

		double GetGradient(ulong i, ulong iVar, FVUtil::COORDINATES dir)
		{
			return GradientPrimMat(i, iVar * (ulong)nDim + (ulong)dir);
		}

		double* GetGradientPrim()
		{
			return GradientPrimMat.block_ptr;
		}
		double* GetGradientPrim(ulong i)
		{
			return GradientPrimMat(i);
		}
		double GetGradientPrim(ulong i, ulong iVar)
		{
			return GradientPrimMat(i, iVar);
		}

		double GetGradientPrim(ulong i, ulong iVar, FVUtil::COORDINATES dir)
		{
			return GradientPrimMat(i, iVar * (ulong)nDim + (ulong)dir);
		}

		double* GetSource()
		{
			return SourceMat.block_ptr;
		}
		double* GetSource(ulong i)
		{
			return SourceMat(i);
		}
		double GetSource(ulong i, ulong iVar)
		{
			return SourceMat(i, iVar);
		}

		double* GetSlopeLimiter()
		{
			return SlopeLimiterMat.block_ptr;
		}
		double* GetSlopeLimiter(ulong i)
		{
			return SlopeLimiterMat(i);
		}
		double GetSlopeLimiter(ulong i, ulong iVar)
		{
			return SlopeLimiterMat(i, iVar);
		}


	};
}