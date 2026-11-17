#pragma once
#include "FluxNumerics.hpp"

namespace FVNumerics
{
	class CRusanovFlux : public CConvectiveFlux
	{
	protected:

	public:
		CRusanovFlux(FVModels::CFluidModel* gas_model, FVModels::CTConductivityModel* therm_model, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
			: CConvectiveFlux(gas_model ,therm_model, mesh, variables)
		{ }
		FVUtil::Vec4 computeFluxAt(ulong edge) override;
	};

	FVUtil::Vec4 CRusanovFlux::computeFluxAt(ulong edge)
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

		//DAVIS WAVESPEED
		const double u_L_abs = std::abs(LState.rho_u() / LState.rho()), u_R_abs = std::abs(RState.rho_u() / RState.rho());
		SetGasModelFromCons(LState);
		const double a_L = GasModel->GetSpeedSound();
		SetGasModelFromCons(RState);
		const double a_R = GasModel->GetSpeedSound();

		const double WaveSpeed = std::max(u_L_abs + a_L, u_R_abs + a_R); //max(|u_L| + a_L, |u_R| + a_R)

		const FluxCons LFlux = GetFluxFromCons(LState), RFlux = GetFluxFromCons(RState);
		flux = 0.5 * (LFlux + RFlux - WaveSpeed * (RState - LState));
		TransformInvT(curr_edge, flux);

		return flux;
	}

}