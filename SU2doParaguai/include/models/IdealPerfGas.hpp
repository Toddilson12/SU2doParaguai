#pragma once

#include "FluidModel.hpp"
#include "../config/config.hpp"
#include "../common.hpp"


namespace FVModels
{
	/*
	Ideal Gas Law + Calorifically Perfect Gas
	*/
	class CIdealPerfGas : CFluidModel
	{
	protected:
		double Gamma;
		double GammaMinOne;
		double GasConstant; // R, not universal

	public:
		CIdealPerfGas(FVSolver::CConfig* config, double R, double gamma);
		void SetThermoState_p_rho(double p, double rho) override;
		void SetThermoState_p_T(double p, double T) override;
		void SetThermoState_h_s(double h, double s) override;
		void SetThermoState_rho_T(double rho, double T) override;
		void SetThermoState_rho_e(double rho, double e) override;
	};

}
