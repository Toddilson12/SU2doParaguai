#pragma once



#include "../models/FluidModel.hpp"
#include "../models/TConductivityModel.hpp"
#include "../mesh/mesh.hpp"
#include "../variables/Variables.hpp"


namespace FVNumerics
{
	/*
	encodePrimitives
	decodePrimitives
	stateVector
	encodeFlux -> depends on models!

	*/


	struct StateCons : FVUtil::Vec4
	{
		double& rho() { return v[0]; }
		double& E() { return v[1]; } // E = e + 0.5 * V^2
		//we could use parameter for diff dimensions
		double& rho_u() { return v[2]; }
		double& rho_v() { return v[3]; }

		double& rho_u_i(ulong i) { return v[2 + i]; }

		//const methods
		const double& rho() const { return v[0]; }
		const double& E() const { return v[1]; }
		const double& rho_u() const { return v[2]; }
		const double& rho_v() const { return v[3]; }
		const double& rho_u_i(ulong i) const { return v[2 + i]; }
	};

	struct StatePrim : FVUtil::Vec4
	{
		double& rho() { return v[0]; }
		double& p() { return v[1]; }
		double& U() { return v[2]; }
		double& V() { return v[3]; }

		double& u_i(ulong i) { return v[2 + i]; } // i = 0 -> x; i = 1 -> y; etc

		//const methods
		const double& rho() const { return v[0]; }
		const double& p() const { return v[1]; }
		const double& U() const { return v[2]; }
		const double& V() const { return v[3]; }

		const double& u_i(ulong i) const { return v[2 + i]; }
	};


	struct FluxCons : FVUtil::Vec4 //2D only
	{
		double& rho_u() { return v[0]; }
		double& En_u() { return v[1]; } 
		double& mom_u() { return v[2]; }
		double& mom_v() { return v[3]; }
		
		double& mom_u_i(ulong i) { return v[2 + i]; }
		//const methods
		const double& rho_u() const { return v[0]; }
		const double& En_u() const { return v[1]; } 
		const double& mom_u() const { return v[2]; }
		const double& mom_v() const { return v[3]; }
		
		const double& mom_u_i(ulong i) const { return v[2 + i]; }
	};

	//Convective Flux
	class CConvectiveFlux
	{
	protected:
		FVModels::CFluidModel* GasModel;
		FVModels::CTConductivityModel* ThermCondModel;
		FVMesh::CMesh* Mesh;
		FVVariables::CVariableContainer* Variables;

		StateCons LState;
		StateCons RState;
		//double NormalAngle; //Angle with x-direction [rad]
		//double FaceArea; // Flux surface area
		void SetGasModelFromCons(const StateCons& state)
		{
			const double u = state.rho_u() / state.rho(), v = state.rho_v() / state.rho();
			const double e_int = state.E() - 0.5 * (u * u + v * v);
			GasModel->SetThermoState_rho_e(state.rho(), e_int);
		}

		void TransformT(const FVMesh::Edge& edge, FVUtil::Vec4& vector) // transforms vector to T(theta) * vector
		{
			const FVUtil::Vec4 temp_vec = vector;
			const double invMag = 1.0 / FVUtil::Vec2::Mag(edge.Normal);
			const double cosTheta = invMag * edge.Normal.x, sinTheta = invMag * edge.Normal.y;

			vector.v[2] = temp_vec.v[2] * cosTheta + temp_vec.v[3] * sinTheta;
			vector.v[3] = temp_vec.v[2] * (-sinTheta) + temp_vec.v[3] * cosTheta;
		}
		void TransformInvT(const FVMesh::Edge& edge, FVUtil::Vec4& vector) // transforms vector to invT(theta) * vector
		{
			const FVUtil::Vec4 temp_vec = vector;
			const double invMag = 1.0 / FVUtil::Vec2::Mag(edge.Normal);
			const double cosTheta = invMag * edge.Normal.x, sinTheta = invMag * edge.Normal.y;

			vector.v[2] = temp_vec.v[2] * cosTheta + temp_vec.v[3] * (-sinTheta);
			vector.v[3] = temp_vec.v[2] * sinTheta + temp_vec.v[3] * cosTheta;
		}

		FluxCons GetFluxFromCons(const StateCons& state)
		{
			FluxCons flux;
			//Set thermo state
			const double u = state.rho_u() / state.rho(), v = state.rho_v() / state.rho();
			const double e_int = state.E() - 0.5 * (u * u + v * v);
			GasModel->SetThermoState_rho_e(state.rho(), e_int);
			const double p = GasModel->GetPressure();

			flux.rho_u() = state.rho_u();
			flux.mom_u() = state.rho_u() * u + p;
			flux.mom_v() = state.rho_u() * v;
			flux.En_u() = u * (state.E() + p);
			
			return flux;
		}

		StatePrim GetPrimitive(const StateCons& state)
		{
			StatePrim prim;
			
			const double u = state.rho_u() / state.rho(), v = state.rho_v() / state.rho();
			const double e_int = state.E() - 0.5 * (u * u + v * v);
			GasModel->SetThermoState_rho_e(state.rho(), e_int);
			const double p = GasModel->GetPressure();
			
			prim.rho() = state.rho();
			prim.p() = p;
			prim.U() = u;
			prim.V() = v;

			return prim;
		}
	
	public:
		CConvectiveFlux(FVModels::CFluidModel* gas_model, FVModels::CTConductivityModel* therm_model, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables);
		virtual FVUtil::Vec4 computeFluxAt(ulong edge) = 0; //computes flux in mesh direction for edge
	};

		


}