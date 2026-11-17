#pragma once
#include "../common.hpp"
#include "../config/config.hpp"
#include "../mesh/mesh.hpp"
#include "../variables/Variables.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>


namespace FVOutput
{
	enum class OutputLabel
	{
		ITERATION,
		TIME
	};
	enum class MeshType
	{
		PRIMAL,
		DUAL
	};


	class COutput
	{
	protected:
		std::string OutputPath;
		std::string	FileNamePrefix;
		OutputLabel LabelSpecific;
		//not used now
		FVSolver::CConfig* Config;

		std::string CurrLabel;

		std::ofstream FileWriter;
	public:
		COutput(FVSolver::CConfig* config, const std::string& output_path, const std::string& filename_prefix, OutputLabel label);
		
		
		template<typename LabelType>
		void SetLabel(LabelType label)
		{
			CurrLabel = std::to_string(label);
		}

		template<MeshType type>
		bool OutputCurrState(FVMesh::CMesh* mesh, FVVariables::CVariableContainer* data, const FVVariables::CIndices* indices);

	};
	
}