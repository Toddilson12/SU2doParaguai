#pragma once

#include "FluxNumerics.hpp"

namespace FVNumerics
{
	//TODO: tranform 1D to 2D
	class CHLLCFlux : public CConvectiveFlux
	{
	protected:
		FluxCons computeFluxStar(const StateCons& U_K, const StatePrim& W_K, double S_K, double S_star)
		{
			/*
			Vec3 U_star_K =
				U_K.rho() * ((S_K - W_K.u()) / (S_K - S_star)) *
				Vec3({
				1.0,
				S_star,
				U_K.E() / U_K.rho() + (S_star - W_K.u()) * (S_star + W_K.p() / (U_K.rho() * (S_K - W_K.u())))
					});
			const FluxFCons F_K = encodeFluxF_ideal_ideal(U_K);
			return F_K + S_K * (U_star_K - U_K);
			*/
			return FluxCons();
		}
	public:
		CHLLCFlux(FVModels::CFluidModel* gas_model, FVModels::CTConductivityModel* therm_model, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
			: CConvectiveFlux(gas_model, therm_model, mesh, variables)
		{
		}

		FVUtil::Vec4 computeFluxAt(ulong edge) override
		{
			FVUtil::Vec4 flux;
			const FVMesh::Edge& curr_edge = Mesh->GetEdge(edge);
			const ulong LNode = curr_edge.nBeg, RNode = curr_edge.nEnd;

			//Get LState and RState (need standard for the whole program)
			//Need coordinate system transform
			LState = StateCons{ {
				Variables->GetSolution(LNode, Variables->Indices.Cons_rho()),
				Variables->GetSolution(LNode, Variables->Indices.Cons_E()),
				Variables->GetSolution(LNode, Variables->Indices.Cons_rhoU()),
				Variables->GetSolution(LNode, Variables->Indices.Cons_rhoV()),
				} };
			RState = StateCons{ {
				Variables->GetSolution(RNode, Variables->Indices.Cons_rho()),
				Variables->GetSolution(RNode, Variables->Indices.Cons_E()),
				Variables->GetSolution(RNode, Variables->Indices.Cons_rhoU()),
				Variables->GetSolution(RNode, Variables->Indices.Cons_rhoV()),
				} };

			TransformT(curr_edge, LState);
			TransformT(curr_edge, RState);

			StatePrim W_L = GetPrimitive(LState);
			StatePrim W_R = GetPrimitive(RState);
			const double gamma = GasModel->GetCp() / GasModel->GetCv();

			const double sqrt_rho_L = std::sqrt(LState.rho());
			const double sqrt_rho_R = std::sqrt(RState.rho());

			const double u_til = (sqrt_rho_L * W_L.U() + sqrt_rho_R * W_R.U()) / (sqrt_rho_L + sqrt_rho_R);

			const double H_L = (LState.E() + W_L.p()) / LState.rho();
			const double H_R = (RState.E() + W_R.p()) / RState.rho();

			const double H_til = (sqrt_rho_L * H_L + sqrt_rho_R * H_R) / (sqrt_rho_L + sqrt_rho_R);
			const double a_til = std::sqrt((gamma - 1.0) * (H_til - 0.5 * u_til * u_til));

			const double S_L = u_til - a_til;
			const double S_R = u_til + a_til;

			const double S_star = (W_R.p() - W_L.p() + LState.rho_u() * (S_L - W_L.U()) - RState.rho_u() * (S_R - W_R.U()))
				/
				(LState.rho() * (S_L - W_L.U()) - RState.rho() * (S_R - W_R.U()));


			// HLLC Flux (Rankie-Hugoniot)

			//SAMPLING
			//First check the simple states (no star region calculations)
			if (0.0 <= S_L)
				flux = GetFluxFromCons(LState);
			else if (0.0 >= S_R)
				flux = GetFluxFromCons(RState);

			//Flux in star region
			else if (S_L < 0.0 && S_star >= 0.0)
				flux = computeFluxStar(LState, W_L, S_L, S_star);
			else if (S_star < 0.0 && S_R >= 0.0)
				flux = computeFluxStar(RState, W_R, S_R, S_star);

			TransformInvT(curr_edge, flux);
			return flux;
		}

	};


}