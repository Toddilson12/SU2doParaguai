#include "../../include/numerics/FluxNumerics.hpp"

namespace FVNumerics
{
	CConvectiveFlux::CConvectiveFlux(FVModels::CFluidModel* gas_model, FVModels::CTConductivityModel* therm_model, FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables)
		: GasModel(gas_model), ThermCondModel(therm_model), Mesh(Mesh), Variables(variables)
	{

	}



}