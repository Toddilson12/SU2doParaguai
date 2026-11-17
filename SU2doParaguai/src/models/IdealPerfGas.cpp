#include "../../include/models/IdealPerfGas.hpp"
#include <cmath>

namespace FVModels
{

	CIdealPerfGas::CIdealPerfGas(FVSolver::CConfig* config, double R, double gamma)
		: CFluidModel(config)
	{
		this->GasConstant = R;
		
		this->Gamma = gamma;
		this->GammaMinOne = Gamma - 1.0;
		this->Cv = GasConstant / (GammaMinOne);
		this->Cp = Cv + GasConstant;
	}


	void CIdealPerfGas::SetThermoState_p_rho(double p, double rho)
	{
		Pressure = p;
		Density = rho;
		
		Temperature = Pressure / (GasConstant * Density);
		SpeedSound = std::sqrt(Gamma * GasConstant * Temperature);
		InternalEnergy = Cv * Temperature;
		Enthalpy = Cp * Temperature;
	}

	void CIdealPerfGas::SetThermoState_p_T(double p, double T)
	{
		Pressure = p;
		Temperature = T;

		Density = Pressure / (GasConstant * Temperature);
		SpeedSound = std::sqrt(Gamma * GasConstant * Temperature);
		InternalEnergy = Cv * Temperature;
		Enthalpy = Cp * Temperature;
	}

	//entropy not used for now
	void CIdealPerfGas::SetThermoState_h_s(double h, double s)
	{
	}

	void CIdealPerfGas::SetThermoState_rho_T(double rho, double T)
	{
		Density = rho;
		Temperature = T;
		
		Pressure = Density * GasConstant * Temperature;
		SpeedSound = std::sqrt(Gamma * GasConstant * Temperature);
		InternalEnergy = Cv * Temperature;
		Enthalpy = Cp * Temperature;
	}
	void CIdealPerfGas::SetThermoState_rho_e(double rho, double e)
	{
		Density = rho;
		InternalEnergy = e;

		Temperature = InternalEnergy / Cv;
		Pressure = Density * GasConstant * Temperature;
		SpeedSound = std::sqrt(Gamma * GasConstant * Temperature);
		Enthalpy = Cp * Temperature;
	}

}