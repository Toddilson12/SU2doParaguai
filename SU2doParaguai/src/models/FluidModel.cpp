#include "../../include/common.hpp"
#include "../../include/models/FluidModel.hpp"


namespace FVModels
{
	CFluidModel::CFluidModel(FVSolver::CConfig* config)
	{
		this->SetViscosityModel(config);
		this->SetTConductivityModel(config);
	}

	CFluidModel::~CFluidModel()
	{
		delete this->ViscosityModel;
		delete this->TConductivityModel;
	}
	
	CViscosityModel* CFluidModel::SetViscosityModel(FVSolver::CConfig* config)
	{
		return nullptr;
	}
	
	CTConductivityModel* CFluidModel::SetTConductivityModel(FVSolver::CConfig* config)
	{
		return nullptr;
	}
}
