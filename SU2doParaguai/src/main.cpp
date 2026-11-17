#include <iostream>
#include <string>

#include "../include/mesh/mesh.hpp"
#include "../include/mesh/meshParser.hpp"
#include "../include/variables/EulerVars.hpp"
#include "../include/output/Output.hpp"
#include "../include/gradients/Gradient.hpp"
#include "../include/gradients/GreenGaussGradient.hpp"
#include "../include/gradients/LeastSquaresGradient.hpp"

#include "../include/numerics/FluxNumerics.hpp"
#include "../include/numerics/RusanovFlux.hpp"
#include "../include/numerics/HLLCFlux.hpp"

double init_cond(double x, double y)
{
	return (3.0 + 2.0 * std::cos(x) + 5.0 * y);
}

int main(int argc, char* argv[])
{
	std::string mesh_path = (argc > 1) ? std::string(argv[1]) : "\0";
	if (mesh_path.back() != '/' && mesh_path.back() != '\\')
		mesh_path += '/';

	std::string mesh_filename = "rect.su2";
	

	FVMesh::CMeshParser mesh_parser(mesh_path + mesh_filename);
	FVMesh::CMesh mesh = mesh_parser.ReadMeshFile();
	mesh.InitDualMesh();

	FVSolver::CConfig config;
	config.OutputCons.push_back("DensityCons");
	config.OutputPrims.push_back("Density");
	config.OutputPrims.push_back("Velocity");
	config.OutputGrad.push_back("DensityCons");
	config.OutputGrad.push_back("E");

	FVVariables::CEulerVariable fields(&config, mesh.GetNumPoints());
	for (ulong i = 0; i < mesh.GetNumPoints(); i++)
	{
		FVUtil::Vec2 pos = mesh.GetNode(i);

		fields.SetSolution(i, fields.Indices.Cons_rho(), init_cond(pos.x, pos.y));
		fields.SetPrimitive(i, fields.Indices.Prim_rho(), init_cond(pos.x, pos.y));
		fields.SetPrimitive(i, fields.Indices.Prim_u(), 1.0);
		fields.SetPrimitive(i, fields.Indices.Prim_v(), 0.0);
	}
	
	//FVGradient::CGreenGaussGradient GreenGaussGrad(&mesh, &fields);
	//GreenGaussGrad.ComputeGradients({ true, true }, 0, 0);
	FVGradient::CLeastSquaresGradient WLS({ true, true }, &mesh, &fields);
	WLS.ComputeGradients({ true, true }, 0, 0);

	FVOutput::COutput output(&config, mesh_path, std::string("test_"), FVOutput::OutputLabel::ITERATION);
	output.SetLabel<int>(0);
	std::cout << "a\n";

	output.OutputCurrState<FVOutput::MeshType::PRIMAL>(&mesh, &fields, &fields.Indices);



	std::cin.get();

	return 0;
}