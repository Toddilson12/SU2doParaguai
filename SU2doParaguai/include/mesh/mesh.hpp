#pragma once
#include "common.hpp"

#include <string>
#include <vector>
#include <map>


namespace FVMesh
{
	enum class ElemType : unsigned int //VTK
	{
		LINE = 3,
		TRIANGLE = 5,
		QUADRILATERAL = 9
	};


	struct Elem
	{
		ElemType type;
		std::vector<int> nodes;
	};

	struct Marker
	{
		std::string	tag;
		std::vector<Elem> elems;
	};


	struct Edge
	{
		int nBeg, nEnd;
		int eLef = -1, eRig = -1; //elements
		FVUtil::Vec2 Normal; // aligned w/ beg -> end
		bool IsBoundary = false;


		bool operator==(const Edge& other) const
		{
			if ((nBeg == other.nBeg && nEnd == other.nEnd) || (nBeg == other.nEnd && nEnd == other.nBeg))
				return true;
			return false;
		}
	};

	struct BoundaryEdge
	{
		int nBeg, nEnd;
		FVUtil::Vec2 Normal; //outward!
	};

	class CMesh
	{
	private:
		int nDime;
		unsigned long nElem, nPoin, nMark;
		//TODO: free unused memory after initialization of mesh
		
		//Original Mesh from parser
		std::vector<FVUtil::Vec2> Nodes;
		std::vector<Elem> Elems;
		std::vector<Marker> Marks;

		//Dual mesh / Prim mesh data
		std::vector<Edge> Edges;
		std::vector<FVUtil::Vec2> Centroids; //idx by elem
		std::vector<double> Volumes; // idx by node
		ulong nIntEdges; // number of internal edges
		ulong nBoundEdges; //number of boundary edges
		std::vector<BoundaryEdge> BoundaryEdges;

		bool isStatic = true;
	public:
		CMesh();
		CMesh(int dime, unsigned long elem, unsigned long poin, unsigned long mark, std::vector<FVUtil::Vec2>&& nodes, std::vector<Elem>&& elems, std::vector<Marker>&& marks, bool isStatic_ = true);

		void InitDualMesh();
		
		bool isMeshStatic() { return isStatic; }

		int GetDimensions() { return nDime; }
		ulong GetNumElem() { return nElem; }
		ulong GetNumPoints() { return nPoin; }
		ulong GetNumMarkers() { return nMark; }
		ulong GetNumEdges() { return Edges.size(); }
		ulong GetNumInternalEdges() { return nIntEdges; }
		ulong GetNumBoundaryEdges() { return nBoundEdges; }
		
		const FVUtil::Vec2* GetNodes_raw() const { return Nodes.data(); }
		const Elem* GetElements_raw() const { return Elems.data(); }
		const Marker* GetMarkers_raw() const { return Marks.data(); }
		const Edge* GetEdges_raw() const { return Edges.data(); }
		const FVUtil::Vec2* GetCentroids_raw() const { return Centroids.data(); }
		const double* GetVolumes_raw() const { return Volumes.data(); }
		const double GetVolumeAt(ulong i) const { return Volumes[i]; }
		const Edge& GetEdge(ulong e) const { return Edges[e]; }
		const BoundaryEdge& GetBoundaryEdge(ulong e) const { return BoundaryEdges[e]; }
		const FVUtil::Vec2& GetNode(ulong i) const { return Nodes[i]; }

	};

}