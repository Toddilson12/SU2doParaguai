#pragma once
#include "../common.hpp"
#include <vector>
#include <string>

namespace FVSolver
{
	class CConfig
	{
	public:
		unsigned short 
			EOSModel,
			ThermalModel,
			InvFluxScheme;
		bool
			MUSCL;
		unsigned short
			SlopeLimiter;
		double
			MachInf,
			ReInf,
			CharLen,
			AoA,
			PressInf,
			DensInf,
			TempInf;
		double
			CFL;
		std::vector<std::string> OutputPrims;
		std::vector<std::string> OutputCons;
		std::vector<std::string> OutputRes;
		std::vector<std::string> OutputGrad;
		std::vector<std::string> OutputGradPrim;
		std::vector<std::string> OutputSource;
		std::vector<std::string> OutputSlopeLimiter;
	};
}