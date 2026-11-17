#pragma once

#include "TConductivityModel.hpp"
#include "ViscosityModel.hpp"
#include "../config/config.hpp"
#include "../common.hpp"

namespace FVModels
{
	class CFluidModel
	{
	protected:
		double
			InternalEnergy,
			Enthalpy,
			SpeedSound,
			Pressure,
			Temperature,
			Density,
			Cp,
			Cv,
			Kt, //thermal conductivity
			Mu; //viscosity
		//unique pointers to visc, therm cond, diff models
		//use static methods for setting visc, therm cond, diff, ... models
		CViscosityModel* ViscosityModel;
		CTConductivityModel* TConductivityModel;

		static CViscosityModel* SetViscosityModel(FVSolver::CConfig* config); 
		static CTConductivityModel* SetTConductivityModel(FVSolver::CConfig* config);
			
	public:
		virtual ~CFluidModel(); //DONT FORGET TO DELETE POINTERS
		CFluidModel(FVSolver::CConfig* config); //generic constructor to grab the pointers' instances

		double GetInternalEnergy() const { return InternalEnergy; }
		double GetEnthalpy() const { return Enthalpy; }
		double GetSpeedSound() const { return SpeedSound; }
		double GetPressure() const { return Pressure; }
		double GetTemperature() const { return Temperature; }
		double GetDensity() const { return Density; }
		double GetCp() const { return Cp; }
		double GetCv() const { return Cv; }
		
		// only calls the pointers' methods
		double GetMu() const // only laminar for now
		{
			//calc stuff and assign to mu
			return Mu;
		}

		double GetKt() const // only laminar for now
		{
			//calc stuff and assign to kt
			return Kt;
		}

		//get lin visc by calc using curr state and funcs from visc model and return mu
		// ...
			
		// Setters (for specific thermostates)
		virtual void SetThermoState_p_rho(double p, double rho) = 0;
		virtual void SetThermoState_p_T(double p, double T) = 0;
		virtual void SetThermoState_h_s(double h, double s) = 0;
		virtual void SetThermoState_rho_T(double rho, double T) = 0;
		virtual void SetThermoState_rho_e(double rho, double e) = 0;
	};
	


}