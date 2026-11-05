#pragma once

#include <cstddef>
#include <stdlib.h>
#include <vector>

//////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Max-flow / min-cut algorithm based on two search trees and orphan adoption
// 
// This code is based on ideas first introduced in:
// Yuri Boykov and Vladimir Kolmogorov "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for 
// Energy Minimization in Vision". IEEE Trans on Pattern Analysis and Machine Intelligence, 26, 2004.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////

#define GC_UNLINKED (-1)
#define GC_TERMINAL (-2)
#define GC_ORPHAN   (-3)

struct GraphNode {
	std::vector<int> edges;	// list of incident edges (indices in the edges array)
	double		terminCap;	// capacity of the terminal edge (positive if connected to source, negative if connected to sink)
	int			parent;		// parent edge in the tree (index in the edges array, or GC_TERMINAL, GC_ORPHAN, GC_UNLINKED)
	int			pdepth;	    // node processing depth (for breadth-first search)
	int			distance;	// distance from the terminal
	bool		isSink;		// indicator showing whether the node is in the source or sink tree (if parent != GC_UNLINKED)
	bool		isActive;	// indicator showing whether the node is in the active queue

	GraphNode() : parent(GC_UNLINKED), pdepth(0), distance(0), isSink(false), terminCap(0), isActive(false) {}

	void SetDistance(int currentDepth, int nodeDistance) {
		pdepth = currentDepth;
		distance = nodeDistance;
	}

	void SetTree(int parentEdge, bool sink, int currentDepth, int nodeDistance) {
		parent = parentEdge;
		isSink = sink;
		SetDistance(currentDepth, nodeDistance);
	}
};

struct GraphEdge {
	int nodes[2];				// endpoints (indices in the nodes array)
	double residualCapacity[2];	// residual capacity in each direction

	GraphEdge(int n1, int n2, double capacity) : nodes{ n1, n2 }, residualCapacity{ capacity, capacity } {}

	int OppositeNode(int nodeIndex) const { return nodes[(nodeIndex == nodes[0])]; }

	double GetResidualCapacity(int fromNode) const { return residualCapacity[(nodes[0] != fromNode)]; }
	void IncrementResidualCapacity(int fromNode, double value) { residualCapacity[(nodes[0] != fromNode)] += value; }
};

class GC_Graph {
	std::vector<GraphNode> nodes; // list of nodes (size = number of nodes)
	std::vector<GraphEdge> edges; // list of edges (size = number of edges)
	std::vector<int> orphanNodes; // orphan nodes (node indices) in the source and sink trees
	std::vector<int> activeQueue; // active nodes (node indices) in the source and sink trees
	size_t activeQueuePosition;   // position of the next node to process in the active queue

	double				flow;		// total flow
	int					pdepth;		// current processing depth

public:
	GC_Graph(int maxNodes)
		: nodes(maxNodes), activeQueuePosition(0), flow(0), pdepth(0) {}

	void AddEdge(int i, int j, double capacity) { edges.emplace_back(i, j, capacity); }

	void AddCapacity(int i, double sourceCapacity, double sinkCapacity) {
		double& d = nodes[i].terminCap;
		const double dd = (d > 0) ? (d) : (0);
		const double csn = sinkCapacity - d;
		flow += dd + std::min(sourceCapacity, csn);
		d += sourceCapacity - sinkCapacity;
	}

	double ComputeMaxFlow() {
		int edge, i, fromNode, node = GC_UNLINKED;

		// Initialize incident edges on nodes
		std::vector<int> incidentEdgesCount(nodes.size(), 0);
		for (i = 0; i < (int)edges.size(); i++) {
			incidentEdgesCount[edges[i].nodes[0]]++;
			incidentEdgesCount[edges[i].nodes[1]]++;
		}
		for (i = 0; i < (int)nodes.size(); i++) nodes[i].edges.reserve(incidentEdgesCount[i]);
		for (i = 0; i < (int)edges.size(); i++) {
			nodes[edges[i].nodes[0]].edges.push_back(i);
			nodes[edges[i].nodes[1]].edges.push_back(i);
		}

		// Initialize active node queue
		activeQueue.reserve(nodes.size());
		for (i = 0; i < (int)nodes.size(); i++) if (nodes[i].terminCap != 0) {
				nodes[i].SetTree(GC_TERMINAL, nodes[i].terminCap < 0, 0, 1);
				PushActiveNode(i);
		}

		orphanNodes.reserve(nodes.size());

		// Main loop
		i = node;
		while (validActiveNode(i)) {
			nodes[i].isActive = false;
			if (GrowTree(i, edge, fromNode)) {
				nodes[i].isActive = true;
				AugmentAndAdopt(edge, fromNode);
			}
			else i = GC_UNLINKED;
		}

		return flow;
	}

	int WhichLabel(int i) const { return (nodes[i].parent != GC_UNLINKED && nodes[i].isSink); }

private:
	int OtherEndpoint(int edgeIndex, int nodeIndex) const { return edges[edgeIndex].OppositeNode(nodeIndex); }

	double ResidualCapacity(int edgeIndex, int fromNode) const { return edges[edgeIndex].GetResidualCapacity(fromNode); }

	void IncrementResidualCapacity(int edgeIndex, int fromNode, double value) { edges[edgeIndex].IncrementResidualCapacity(fromNode, value); }

	bool validActiveNode(int& i) {
		if (i == GC_UNLINKED) { if ((i = NextActiveNode()) == GC_UNLINKED) return false; }
		if (nodes[i].parent == GC_UNLINKED) { if ((i = NextActiveNode()) == GC_UNLINKED) return false; }
		return true;
	}

	void PushActiveNode(int i) {
		if (!nodes[i].isActive)
		{
			activeQueue.push_back(i);
			nodes[i].isActive = true;
		}
	}

	int NextActiveNode() {
		while (activeQueuePosition != activeQueue.size()) {
			const int i = activeQueue[activeQueuePosition++];
			nodes[i].isActive = false;
			if (nodes[i].parent != GC_UNLINKED) return i;
		}
		return GC_UNLINKED;
	}

	void PushOrphan(int i) {
		nodes[i].parent = GC_ORPHAN;
		orphanNodes.push_back(i);
	}

	bool GrowTree(int r, int& edge, int& fromNode) {
		const GraphNode& nr = nodes[r];
		const bool ss = nr.isSink;
		edge = fromNode = GC_UNLINKED;

		for (int edgeIndex : nr.edges) {
			const int i = OtherEndpoint(edgeIndex, r);
			int ir = (ss) ? (i) : (r);
			if (ResidualCapacity(edgeIndex, ir) != 0) {
				GraphNode& ni = nodes[i];
				if (ni.parent == GC_UNLINKED) {
					ni.SetTree(edgeIndex, ss, nr.pdepth, nr.distance + 1);
					PushActiveNode(i);
				}
				else if (ni.isSink ^ ss) { edge = edgeIndex; fromNode = ir; break; }
				else if (ni.pdepth <= nr.pdepth && ni.distance > nr.distance) {
					ni.parent = edgeIndex;
					ni.SetDistance(nr.pdepth, nr.distance + 1);
				}
			}
		}
		pdepth++;

		return (edge != GC_UNLINKED);
	}

	void AugmentAndAdopt(int centralEdge, int fromNode) {
		const int toNode = OtherEndpoint(centralEdge, fromNode);
		double capacity, min_cap = ResidualCapacity(centralEdge, fromNode);
		int i, edge, parent;

		for (i = fromNode; (edge = nodes[i].parent) != GC_TERMINAL; ) {
			i = OtherEndpoint(edge, i);
			capacity = ResidualCapacity(edge, i);
			if (min_cap > capacity) min_cap = capacity;
		}
		if (min_cap > nodes[i].terminCap) min_cap = nodes[i].terminCap;

		for (i = toNode; (edge = nodes[i].parent) != GC_TERMINAL;) {
			capacity = ResidualCapacity(edge, i);
			i = OtherEndpoint(edge, i);
			if (min_cap > capacity) min_cap = capacity;
		}
		if (min_cap > -nodes[i].terminCap) min_cap = -nodes[i].terminCap;

		IncrementResidualCapacity(centralEdge, toNode, min_cap);
		IncrementResidualCapacity(centralEdge, fromNode, -min_cap);
		for (i = fromNode; (edge = nodes[i].parent) != GC_TERMINAL; i = parent) {
			parent = OtherEndpoint(edge, i);
			IncrementResidualCapacity(edge, i, min_cap);
			IncrementResidualCapacity(edge, parent, -min_cap);
			if (ResidualCapacity(edge, parent) == 0) PushOrphan(i);
		}
		nodes[i].terminCap -= min_cap;
		if (nodes[i].terminCap == 0) PushOrphan(i);

		for (i = toNode; (edge = nodes[i].parent) != GC_TERMINAL; i = parent) {
			parent = OtherEndpoint(edge, i);
			IncrementResidualCapacity(edge, parent, min_cap);
			IncrementResidualCapacity(edge, i, -min_cap);
			if (ResidualCapacity(edge, i) == 0) PushOrphan(i);
		}
		nodes[i].terminCap += min_cap;
		if (nodes[i].terminCap == 0) PushOrphan(i);

		for (i = 0; i < (int)orphanNodes.size(); i++) AdoptOrphan(orphanNodes[i]);
		orphanNodes.clear();

		flow += min_cap;
	}

	void AdoptOrphan(int n) {
		int i, edge, dist, minDist = INT_MAX, minEdge = GC_UNLINKED;
		const GraphNode& ni = nodes[n];
		const bool source = !ni.isSink;

		for (const int edgeIndex : ni.edges) {
			const int g = OtherEndpoint(edgeIndex, n);
			i = (source) ? (g) : (n);
			if (ResidualCapacity(edgeIndex, i) != 0 && (nodes[g].isSink ^ source) && nodes[g].parent != GC_UNLINKED) {
				for (i = g, dist = 0;; dist++, i = OtherEndpoint(edge, i)) {
					if (nodes[i].pdepth == pdepth) { dist += nodes[i].distance; break; }
					edge = nodes[i].parent;
					if (edge == GC_TERMINAL) { nodes[i].SetDistance(pdepth, 1); break; }
					if (edge == GC_ORPHAN) { dist = INT_MAX; break; }
				}
				if (dist < INT_MAX) {
					if (dist < minDist) { minEdge = edgeIndex; minDist = dist; }
					for (i = g; nodes[i].pdepth != pdepth; i = OtherEndpoint(nodes[i].parent, i))
						nodes[i].SetDistance(pdepth, dist--);
				}
			}
		}

		nodes[n].parent = minEdge;
		if (minEdge != GC_UNLINKED) nodes[n].SetDistance(pdepth, minDist + 1);
		else {
			for (const int edgeIndex : ni.edges) {
				i = OtherEndpoint(edgeIndex, n);
				edge = nodes[i].parent;
				if ((nodes[i].isSink ^ source) && edge != GC_UNLINKED) {
					if (ResidualCapacity(edgeIndex, n)) PushActiveNode(i);
					if (edge != GC_TERMINAL && edge != GC_ORPHAN && OtherEndpoint(edge, i) == n) PushOrphan(i);
				}
			}
		}
	}
};
