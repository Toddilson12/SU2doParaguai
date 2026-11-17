#pragma once
#include "mesh.hpp"
#include <string>

namespace FVMesh
{
	class CMeshParser
	{
	private:
		std::string MeshPath;
	public:
		CMeshParser(const std::string& file_path = "");
		CMesh ReadMeshFile(const std::string& file_path = "");
	};
}