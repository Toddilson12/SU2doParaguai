#include "GreenGaussGradient.hpp"



namespace FVGradient
{
	void CGreenGaussGradient::ComputeGradients(const GradOptions& opts, ulong iVarBegin, ulong iVarEnd)
	{
		ulong nPoin = Mesh->GetNumPoints();
		ulong nEdge = Mesh->GetNumEdges();
		ulong nBoundEdge = Mesh->GetNumBoundaryEdges();
		ulong nVars = Variables->GetNumVars();
		ulong nPrimsGrad = Variables->GetNumPrimsGrad();

		if (iVarEnd == -1)
			iVarEnd = nVars - 1;

		// Clear gradient first
		for (ulong i = 0; i < nPoin; i++)
		{
			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
				Variables->SetGradient(i, iVar, { 0.0, 0.0 });
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nPrimsGrad; iPrim++)
					Variables->SetGradientPrim(i, iPrim, { 0.0, 0.0 });
			}
		}

		for (ulong e = 0; e < nEdge; e++)
		{
			const FVMesh::Edge currEdge = Mesh->GetEdge(e);
			const ulong iBeg = currEdge.nBeg, iEnd = currEdge.nEnd;

			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
			{
				FVUtil::Vec2 sum = currEdge.Normal * (Variables->GetSolution(iBeg, iVar) + Variables->GetSolution(iEnd, iVar));
				Variables->AddGradient(iBeg, iVar, sum);
				Variables->AddGradient(iEnd, iVar, sum * (-1.0));
			}
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nPrimsGrad; iPrim++)
				{
					FVUtil::Vec2 sum = currEdge.Normal * (Variables->GetPrimitive(iBeg, iPrim) + Variables->GetPrimitive(iEnd, iPrim));
					Variables->AddGradientPrim(iBeg, iPrim, sum);
					Variables->AddGradientPrim(iEnd, iPrim, sum * (-1.0));
				}
			}
		}

		// Boundary edges
		for (ulong e = 0; e < nBoundEdge; e++)
		{
			const FVMesh::BoundaryEdge currEdge = Mesh->GetBoundaryEdge(e);
			const ulong iBeg = currEdge.nBeg, iEnd = currEdge.nEnd;

			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
			{
				FVUtil::Vec2 sum_ij = currEdge.Normal *
					((5.0 * Variables->GetSolution(iBeg, iVar) + Variables->GetSolution(iEnd, iVar)) / 6.0);

				FVUtil::Vec2 sum_ji = currEdge.Normal *
					((Variables->GetSolution(iBeg, iVar) + 5.0 * Variables->GetSolution(iEnd, iVar)) / 6.0);

				Variables->AddGradient(iBeg, iVar, sum_ij);
				Variables->AddGradient(iEnd, iVar, sum_ji);
			}
			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nPrimsGrad; iPrim++)
				{
					FVUtil::Vec2 sum_ij = currEdge.Normal *
						((5.0 * Variables->GetPrimitive(iBeg, iPrim) + Variables->GetPrimitive(iEnd, iPrim)) / 6.0);

					FVUtil::Vec2 sum_ji = currEdge.Normal *
						((Variables->GetPrimitive(iBeg, iPrim) + 5.0 * Variables->GetPrimitive(iEnd, iPrim)) / 6.0);

					Variables->AddGradientPrim(iBeg, iPrim, sum_ij);
					Variables->AddGradientPrim(iEnd, iPrim, sum_ji);
				}
			}
		}

		for (ulong i = 0; i < nPoin; i++) // multiply by each cell volume
		{
			const double inv_2Vol = 0.5 / Mesh->GetVolumeAt(i);
			for (ulong iVar = iVarBegin; iVar <= iVarEnd; iVar++)
				Variables->ScaleGradient(i, iVar, inv_2Vol);

			if (opts.computePrim)
			{
				for (ulong iPrim = 0; iPrim < nPrimsGrad; iPrim++)
					Variables->ScaleGradientPrim(i, iPrim, inv_2Vol);
			}
		}
	}
}
