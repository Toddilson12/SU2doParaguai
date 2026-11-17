#include "mesh/meshParser.hpp"
#include <sstream>
#include <fstream>
#include <iostream>
#include <utility>

namespace FVMesh
{
	enum class ReaderState
	{
		NDIME,
		NPOIN,
		NELEM,
		NMARK,
		NMARKER_TAG,
		NMARKER_ELEMS,
		NONE
	};

	const std::map<ReaderState, std::string> StateMap =
	{
		{ReaderState::NDIME, "NDIME"},
		{ReaderState::NPOIN, "NPOIN"},
		{ReaderState::NELEM, "NELEM"},
		{ReaderState::NMARK, "NMARK"},
		{ReaderState::NMARKER_TAG, "NMARKER_TAG"},
		{ReaderState::NMARKER_ELEMS, "NMARKER_ELEMS"}
	};

	ReaderState checkLineReaderState(const std::string& line)
	{
		for (const auto& pair : StateMap)
		{
			if (line.find(pair.second) == std::size_t(0))
				return pair.first;
		}
		return ReaderState::NONE;
	}

	int readLineIntVal(const std::string& line)
	{
		size_t eq_pos = line.find('=');
		if (eq_pos != std::string::npos)
			return std::stoi(line.substr(eq_pos + 1));
		return -1;
	}

	std::string readLineStrVal(const std::string& line)
	{
		size_t eq_pos = line.find('=');
		std::string val = "";
		if (eq_pos != std::string::npos)
			val = line.substr(eq_pos + 1);
		
		std::stringstream ss(val);
		std::string	res = "";

		ss >> res;

		return res;
	}





	CMeshParser::CMeshParser(const std::string& file_path)
	{
		MeshPath = file_path;
	}
	
	CMesh CMeshParser::ReadMeshFile(const std::string& file_path)
	{
		MeshPath = (file_path == "") ? MeshPath : file_path;
		
		std::ifstream mesh_file(MeshPath);
		
		std::cout << "Reading mesh file from " << MeshPath << "\n";

		if (!mesh_file.is_open())
		{
			std::cout << "Couldnt open mesh file\n";
			exit(-1);
		}
		

		
		std::string line;

		int nDime = 0;
		int nElem = 0;
		int nPoin = 0;
		int nMark = 0;

		
		std::vector<FVUtil::Vec2> Nodes;
		std::vector<Elem> Elems;
		std::vector<Marker> Marks;
		
		ReaderState rState = ReaderState::NONE;
		while (std::getline(mesh_file, line))
		{
			//Skip empty line
			if (line.empty())
				continue;
			
			if (rState == ReaderState::NONE)
				rState = checkLineReaderState(line);
			
			switch (rState)
			{
			case ReaderState::NDIME:
				nDime = readLineIntVal(line);
				rState = ReaderState::NONE;
				break;
			case ReaderState::NPOIN:
				nPoin = readLineIntVal(line);
				for (size_t n = 0; n < nPoin; n++)
				{
					std::getline(mesh_file, line); 
					std::stringstream ss(line);
					int id; //unused for now, since they should be in order
					double x, y;
					if (ss >> x >> y >> id)
						Nodes.emplace_back(x, y);
				}
				rState = ReaderState::NONE;
				break;
			case ReaderState::NELEM:
				nElem = readLineIntVal(line);
				for (size_t n = 0; n < nElem; n++)
				{
					std::getline(mesh_file, line); 
					std::stringstream ss(line);
					int id; //unused for now, since they should be in order
					double x, y;
					Elem temp;
					if (ss >> id)
					{
						//what a stupid code with lots of memory alloc.
						int val;
						temp.type = (ElemType)id;
						while (ss >> val)
							temp.nodes.push_back(val);
						temp.nodes.pop_back();
						Elems.push_back(temp);
					}
				}

				rState = ReaderState::NONE;
				break;
			case ReaderState::NMARK:
				nMark = readLineIntVal(line);
				Marks.reserve(nMark);

				for (size_t n = 0; n < nMark; n++)
				{
					Marker marker;
					std::getline(mesh_file, line);
					marker.tag = readLineStrVal(line);


					size_t marker_elems = 0;
					std::getline(mesh_file, line);
					marker_elems = readLineIntVal(line);

					for (size_t i = 0; i < marker_elems; i++)
					{
						Elem elem_temp;
						std::getline(mesh_file, line);
						std::stringstream ss(line);
						int id, node;

						if (ss >> id)
						{
							elem_temp.type = (ElemType)id;
							while (ss >> node)
								elem_temp.nodes.push_back(node);
						}

						marker.elems.push_back(elem_temp);
					}

					Marks.push_back(marker);
				}
				rState = ReaderState::NONE;
				break;
			default:
			case ReaderState::NONE:
				break;
			}
			
		}
		
		std::cout << "Done reading mesh file (" << MeshPath << ")\n"
			<< "Mesh data: " << "\n"
			<< "\tDimensions: " << nDime << "\n"
			<< "\tElements: " << nElem << "\n"
			<< "\tPoints/Nodes: " << nPoin << "\n";
		std::cout << "\tMarkers: " << nMark << "\n";
		for (const auto& marker : Marks)
			std::cout << "\t\t\"" << marker.tag << "\": " << marker.elems.size() << " element(s)\n";


		return CMesh(nDime, nElem, nPoin, nMark, std::move(Nodes), std::move(Elems), std::move(Marks));
	}

}