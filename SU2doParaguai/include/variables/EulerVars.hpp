#include "Variables.hpp"
#include "Indices.hpp"
/*
FOR NOW, ONLY 2D IS SUPPORTED. DIMENSIONS ARE HARD CODED TO 2D
*/
namespace FVVariables
{
	class CEulerVariable : public CVariableContainer
	{
	protected:
		//CIndices Indices;
	public:
		CEulerVariable(FVSolver::CConfig* config, ulong n_poin);
		virtual ~CEulerVariable();
		// Setters for the conservative
		
		void SetRho(ulong i, double rho) { SetSolution(i, Indices.Cons_rho(), rho); SetPrimitive(i, Indices.Prim_rho(), rho); }
		void SetRhoField(double rho_0) 
		{
			for (ulong i = 0; i < nPoin; i++)
			{
				SetSolution(i, Indices.Cons_rho(), rho_0);
				SetPrimitive(i, Indices.Prim_rho(), rho_0);
			}
		}
	 
		void SetRhoU(ulong i, double rho_u) { SetSolution(i, Indices.Cons_rhoU(), rho_u); }
		void SetRhoUField(double rho_u_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Cons_rhoU(), rho_u_0);
		}

		void SetRhoV(ulong i, double rho_v) { SetSolution(i, Indices.Cons_rhoV(), rho_v); }
		void SetRhoVField(double rho_v_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Cons_rhoV(), rho_v_0);
		}

		void SetE(ulong i, double E) { SetSolution(i, Indices.Cons_E(), E); }
		void SetEField(double E_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Cons_E(), E_0);	
		}

		// Setters for the primitive
		
		void SetU(ulong i, double u) { SetSolution(i, Indices.Prim_u(), u); }
		void SetUField(double u_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Prim_u(), u_0);
		}

		void SetV(ulong i, double v) { SetSolution(i, Indices.Prim_v(), v); }
		void SetVField(double v_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Prim_v(), v_0);
		}

		void SetP(ulong i, double p) { SetSolution(i, Indices.Prim_p(), p); }
		void SetPField(double p_0)
		{
			for (ulong i = 0; i < nPoin; i++)
				SetSolution(i, Indices.Prim_p(), p_0);	
		}


		// Getters for the conservative
		double GetRho(ulong i) { return GetSolution(i, Indices.Cons_rho()); }
		double GetRhoU(ulong i) { return GetSolution(i, Indices.Cons_rhoU()); }
		double GetRhoV(ulong i) { return GetSolution(i, Indices.Cons_rhoV()); }
		double GetE(ulong i) { return GetSolution(i, Indices.Cons_E()); }

		// Getters for the primitive
		double GetU(ulong i) { return GetPrimitive(i, Indices.Prim_u()); }
		double GetV(ulong i) { return GetPrimitive(i, Indices.Prim_v()); }
		double GetP(ulong i) { return GetPrimitive(i, Indices.Prim_p()); }

		// Getters for additional values (may be create a secondary matrix like SU2 for storing these stuff or add more stuff in the primitives)
		
	};



}
