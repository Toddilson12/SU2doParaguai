#pragma once
#include "../common.hpp"
#include <vector>
#include <string>
#include <unordered_map>

namespace FVVariables
{
	class CIndices
	{
		ulong
			ConsRhoIdx = 0,
			ConsRhoUIdx = 1,
			ConsRhoVIdx = 2,
			ConsEIdx = 3;
		ulong
			PrimRhoIdx = 0,
			PrimUIdx = 1,
			PrimVIdx = 2,
			PrimPIdx = 3;

		std::unordered_map<std::string, ulong> PrimIndexNames;
		std::unordered_map<std::string, ulong> ConsIndexNames;
			
	public:
		CIndices()
		{
			PrimIndexNames["Density"] = 0;
			PrimIndexNames["u"] = 1;
			PrimIndexNames["v"] = 2;
			PrimIndexNames["Pressure"] = 3;

			ConsIndexNames["DensityCons"] = 0;
			ConsIndexNames["rhoU"] = 1;
			ConsIndexNames["rhoV"] = 2;
			ConsIndexNames["E"] = 3;
		}
		//Setters
		void SetConsRhoIdx(ulong idx) { ConsRhoIdx = idx; ConsIndexNames["DensityCons"] = idx; }
		void SetConsRhoUIdx(ulong idx) { ConsRhoUIdx = idx; ConsIndexNames["rhoU"] = idx; }
		void SetConsRhoVIdx(ulong idx) { ConsRhoVIdx = idx; ConsIndexNames["rhoV"] = idx; }
		void SetConsEIdx(ulong idx) { ConsEIdx = idx; ConsIndexNames["E"] = idx; }

		void SetPrimRhoIdx(ulong idx) { PrimRhoIdx = idx; PrimIndexNames["Density"] = idx; }
		void SetPrimUIdx(ulong idx) { PrimUIdx = idx; PrimIndexNames["u"] = idx; }
		void SetPrimVIdx(ulong idx) { PrimVIdx = idx; PrimIndexNames["v"] = idx; }
		void SetPrimPIdx(ulong idx) { PrimPIdx = idx; PrimIndexNames["Pressure"] = idx; }
		//Getters
		ulong Cons_rho() const { return ConsRhoIdx; } 
		ulong Cons_rhoU() const { return ConsRhoUIdx; } 
		ulong Cons_rhoV() const { return ConsRhoVIdx; } 
		ulong Cons_E() const { return ConsEIdx; } 

		ulong Prim_rho() const { return PrimRhoIdx; }
		ulong Prim_u() const { return PrimUIdx; } 
		ulong Prim_v() const { return PrimVIdx; } 
		ulong Prim_p() const { return PrimPIdx; } 
		
		const std::unordered_map<std::string, ulong>& GetPrimIndexNames() const
		{
			return PrimIndexNames;
		}

		const std::unordered_map<std::string, ulong>& GetConsIndexNames() const
		{
			return ConsIndexNames;
		}

		ulong GetPrimIndex(const std::string& name) const
		{
			return PrimIndexNames.at(name);
		}
		ulong GetConsIndex(const std::string& name) const
		{
			return ConsIndexNames.at(name);
		}



	};

};