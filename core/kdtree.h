
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_KDTREE_H
#define PBRT_KDTREE_H
// kdtree.h*
#include "pbrt.h"
#include "geometry.h"
// KdTree Declarations
struct KdNode {
	void init(float p, u_int a) {
		splitPos = p;
		splitAxis = a;
		rightChild = ~(0U);
		hasLeftChild = 0;
	}
	void initLeaf() {
		splitAxis = 3;
		rightChild = ~(0U);
		hasLeftChild = 0;
	}
	// KdNode Data
	float splitPos;
	u_int splitAxis:2;
	u_int hasLeftChild:1;
	u_int rightChild:29;
};
template <class NodeData, class LookupProc> class KdTree {
public:
	// KdTree Public Methods
	KdTree(const vector<NodeData> &data);
	~KdTree() {
		FreeAligned(nodes);
		delete[] nodeData;
	}
	void recursiveBuild(u_int nodeNum, int start, int end,
		vector<const NodeData *> &buildNodes);
	void Lookup(const Point &p, const LookupProc &process,
			float &maxDistSquared) const;
private:
	// KdTree Private Methods
	void privateLookup(u_int nodeNum, const Point &p,
		const LookupProc &process, float &maxDistSquared) const;
	// KdTree Private Data
	KdNode *nodes;
	NodeData *nodeData;
	u_int nNodes, nextFreeNode;
};
template<class NodeData> struct CompareNode {
	CompareNode(int a) { axis = a; }
	int axis;
	bool operator()(const NodeData *d1,
			const NodeData *d2) const {
		return d1->p[axis] == d2->p[axis] ? (d1 < d2) :
			d1->p[axis] < d2->p[axis];
	}
};
// KdTree Method Definitions
template <class NodeData, class LookupProc>
KdTree<NodeData,
       LookupProc>::KdTree(const vector<NodeData> &d) {
	nNodes = d.size();
	nextFreeNode = 1;
	nodes = (KdNode *)AllocAligned(nNodes * sizeof(KdNode));
	nodeData = new NodeData[nNodes];
	vector<const NodeData *> buildNodes;
	for (u_int i = 0; i < nNodes; ++i)
		buildNodes.push_back(&d[i]);
	// Begin the KdTree building process
	recursiveBuild(0, 0, nNodes, buildNodes);
}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::recursiveBuild(u_int nodeNum,
		int start, int end,
		vector<const NodeData *> &buildNodes) {
	// Create leaf node of kd-tree if we've reached the bottom
	if (start + 1 == end) {
		nodes[nodeNum].initLeaf();
		nodeData[nodeNum] = *buildNodes[start];
		return;
	}
	// Choose split direction and partition data
	// Compute bounds of data from _start_ to _end_
	BBox bound;
	for (int i = start; i < end; ++i)
		bound = Union(bound, buildNodes[i]->p);
	int splitAxis = bound.MaximumExtent();
	int splitPos = (start+end)/2;
	std::nth_element(&buildNodes[start], &buildNodes[splitPos],
		&buildNodes[end-1]+1, CompareNode<NodeData>(splitAxis));
	// Allocate kd-tree node and continue recursively
	nodes[nodeNum].init(buildNodes[splitPos]->p[splitAxis],
		splitAxis);
	nodeData[nodeNum] = *buildNodes[splitPos];
	if (start < splitPos) {
		nodes[nodeNum].hasLeftChild = 1;
		u_int childNum = nextFreeNode++;
		recursiveBuild(childNum, start, splitPos, buildNodes);
	}
	if (splitPos+1 < end) {
		nodes[nodeNum].rightChild = nextFreeNode++;
		recursiveBuild(nodes[nodeNum].rightChild, splitPos+1,
		               end, buildNodes);
	}
}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::Lookup(const Point &p,
		const LookupProc &proc,
		float &maxDistSquared) const {
	privateLookup(0, p, proc, maxDistSquared);
}
template <class NodeData, class LookupProc> void
KdTree<NodeData, LookupProc>::privateLookup(u_int nodeNum,
		const Point &p,	const LookupProc &process,
		float &maxDistSquared) const {
	KdNode *node = &nodes[nodeNum];
	// Process kd-tree node's children
	int axis = node->splitAxis;
	if (axis != 3) {
		float dist2 = (p[axis] - node->splitPos) *
			(p[axis] - node->splitPos);
		if (p[axis] <= node->splitPos) {
			if (node->hasLeftChild)
				privateLookup(nodeNum+1, p,
				              process, maxDistSquared);
			if (dist2 < maxDistSquared &&
			    node->rightChild < nNodes)
				privateLookup(node->rightChild,
				              p,
							  process,
							  maxDistSquared);
		}
		else {
			if (node->rightChild < nNodes)
				privateLookup(node->rightChild,
				              p,
							  process,
							  maxDistSquared);
			if (dist2 < maxDistSquared && node->hasLeftChild)
				privateLookup(nodeNum+1,
				              p,
							  process,
							  maxDistSquared);
		}
	}
	// Hand kd-tree node to processing function
	float dist2 = DistanceSquared(nodeData[nodeNum].p, p);
	if (dist2 < maxDistSquared)
		process(nodeData[nodeNum], dist2, maxDistSquared);
}
#endif // PBRT_KDTREE_H
