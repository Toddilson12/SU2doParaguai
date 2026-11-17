#include "../../include/output/Output.hpp"

namespace FVOutput
{
	COutput::COutput(FVSolver::CConfig* config, const std::string& output_path, const std::string& filename_prefix, OutputLabel label)
		: Config(config), OutputPath(output_path), FileNamePrefix(filename_prefix), LabelSpecific(label)
	{
		//correct output path and prefix string
		if (OutputPath.back() != '/' && OutputPath.back() != '\\')
			OutputPath += '/';
		//if (FileNamePrefix.back() != '/' && FileNamePrefix.back() != '\\')
		//	FileNamePrefix += '/';
	}
	template<>
	bool COutput::OutputCurrState<MeshType::PRIMAL>(FVMesh::CMesh* mesh, FVVariables::CVariableContainer* variables, const FVVariables::CIndices* indices)
	{
		const std::string filename = OutputPath + FileNamePrefix + CurrLabel + ".vtu";

		if (FileWriter.is_open())
			FileWriter.close();

		FileWriter.open(filename);

		if (!FileWriter.is_open())
		{
			std::cout << "COutput::OutputCurrState(...): Could not open output file (" << filename << ")\n";
			return false;
		}

		const ulong nPoin = mesh->GetNumPoints();
		const ulong nElem = mesh->GetNumElem();

		//header for xml
		FileWriter << "<?xml version=\"1.0\"?>\n";

		//Start document level element
		FileWriter << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

		//Start UnstructuredGrid element
		FileWriter << "\t<UnstructuredGrid>\n";

		//Start Piece Element
		FileWriter << "\t\t<Piece NumberOfPoints=\"" << nPoin << "\" NumberOfCells=\"" << nElem << "\">\n";

		//start Points
		FileWriter << "\t\t\t<Points>\n";
		FileWriter << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		FileWriter << "\t\t\t\t\t";

		for (ulong i = 0; i < nPoin; i++)
		{
			FileWriter << mesh->GetNodes_raw()[i].x << " "
				<< mesh->GetNodes_raw()[i].y << " 0.0 ";
		}

		FileWriter << "\n";
		FileWriter << "\t\t\t\t</DataArray>\n";
		FileWriter << "\t\t\t</Points>\n";

		//start PointData element
		FileWriter << "\t\t\t<PointData>\n";

		//Solution
		for (const auto& name : Config->OutputCons)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << name << "\"";

			if (name == "RhoVel")
				FileWriter << " NumberOfComponents=\"3\"";

			FileWriter << " format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				if (name == "RhoVel")
				{
					FileWriter << variables->GetSolution(i, indices->Cons_rhoU()) << " "
						<< variables->GetSolution(i, indices->Cons_rhoV()) << " 0.0 ";
				}
				else
				{
					FileWriter << variables->GetSolution(i, indices->GetConsIndex(name)) << " ";
				}
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//Prim
		for (const auto& name : Config->OutputPrims)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << name << "\"";

			if (name == "Velocity")
				FileWriter << " NumberOfComponents=\"3\"";

			FileWriter << " format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				if (name == "Velocity")
				{
					FileWriter << variables->GetPrimitive(i, indices->Prim_u()) << " "
						<< variables->GetPrimitive(i, indices->Prim_v()) << " 0.0 ";
				}
				else
				{
					FileWriter << variables->GetPrimitive(i, indices->GetPrimIndex(name)) << " ";
				}
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//Residual
		for (const auto& name : Config->OutputRes)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Res_" << name << "\" format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				FileWriter << variables->GetResidual(i, indices->GetConsIndex(name)) << " ";
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//Gradient
		for (const auto& name : Config->OutputGrad)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Grad_" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				FileWriter << variables->GetGradient(i, indices->GetConsIndex(name), FVUtil::COORDINATES::X) << " "
					<< variables->GetGradient(i, indices->GetConsIndex(name), FVUtil::COORDINATES::Y) << " 0.0 ";
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//GradientPrim
		for (const auto& name : Config->OutputGradPrim)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Grad_" << name << "\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				FileWriter << variables->GetGradientPrim(i, indices->GetConsIndex(name), FVUtil::COORDINATES::X) << " "
					<< variables->GetGradientPrim(i, indices->GetConsIndex(name), FVUtil::COORDINATES::Y) << " 0.0 ";
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//Source
		for (const auto& name : Config->OutputSource)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Source_" << name << "\" format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				FileWriter << variables->GetSource(i, indices->GetConsIndex(name)) << " ";
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//SlopeLimiter
		for (const auto& name : Config->OutputSlopeLimiter)
		{
			FileWriter << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Slope_" << name << "\" format=\"ascii\">\n\t\t\t\t\t";

			for (ulong i = 0; i < nPoin; i++)
			{
				FileWriter << variables->GetSlopeLimiter(i, indices->GetPrimIndex(name)) << " ";
			}

			FileWriter << "\n\t\t\t\t</DataArray>\n";
		}

		//Finish PointData Element
		FileWriter << "\t\t\t</PointData>\n";

		//Start Cell Element
		FileWriter << "\t\t\t<Cells>\n";

		//Connectivity
		FileWriter << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n\t\t\t\t\t";

		for (ulong el = 0; el < nElem; el++)
		{
			const ulong node_num = mesh->GetElements_raw()[el].nodes.size();

			for (ulong n = 0; n < node_num; n++)
				FileWriter << mesh->GetElements_raw()[el].nodes[n] << " ";
		}

		FileWriter << "\n\t\t\t\t</DataArray>\n";

		//Offsets
		FileWriter << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t";

		ulong cummulative_idx = 0;

		for (ulong el = 0; el < nElem; el++)
		{
			cummulative_idx += mesh->GetElements_raw()[el].nodes.size();
			FileWriter << cummulative_idx << " ";
		}

		FileWriter << "\n\t\t\t\t</DataArray>\n";

		//VTK types
		FileWriter << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t";

		for (ulong el = 0; el < nElem; el++)
		{
			FileWriter << static_cast<unsigned int>(mesh->GetElements_raw()[el].type) << " ";
		}

		FileWriter << "\n\t\t\t\t</DataArray>\n";

		//Finish Cell Element
		FileWriter << "\t\t\t</Cells>\n";

		//Finish Piece Element
		FileWriter << "\t\t</Piece>\n";

		//Finish UnstructuredGrid element
		FileWriter << "\t</UnstructuredGrid>\n";

		//Finish document level element
		FileWriter << "</VTKFile>\n";

		FileWriter.close();

		return true;
	}
}