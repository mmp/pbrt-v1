
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

// kdtree.cpp*
#include "pbrt.h"
#include "primitive.h"
// KdAccelNode Declarations
struct MailboxPrim {
	MailboxPrim(const Reference<Primitive> &p) {
		primitive = p;
		lastMailboxId = -1;
	}
	Reference<Primitive> primitive;
	int lastMailboxId;
};
struct KdAccelNode {
	// KdAccelNode Methods
	void initLeaf(int *primNums, int np,
			MailboxPrim *mailboxPrims, MemoryArena &arena) {
		// Update kd leaf node allocation statistics
		static StatsCounter numLeafMade("Kd-Tree Accelerator",
		                                "Leaf kd-tree nodes made");
		//static StatsCounter maxDepth("Kd-Tree Accelerator",
		//                             "Maximum kd-tree depth");
		static StatsCounter maxLeafPrims("Kd-Tree Accelerator",
			"Maximum number of primitives in leaf node");
		++numLeafMade;
		//maxDepth.Max(depth);
		maxLeafPrims.Max(np);
		static StatsRatio leafPrims("Kd-Tree Accelerator",
			"Avg. number of primitives in leaf nodes");
		leafPrims.Add(np, 1);
		nPrims = np << 2;
		flags |= 3;
		// Store _MailboxPrim *_s for leaf node
		if (np == 0)
			onePrimitive = NULL;
		else if (np == 1)
			onePrimitive = &mailboxPrims[primNums[0]];
		else {
			primitives = (MailboxPrim **)arena.Alloc(np *
				sizeof(MailboxPrim *));
			for (int i = 0; i < np; ++i)
				primitives[i] = &mailboxPrims[primNums[i]];
		}
	}
	void initInterior(int axis, float s) {
		static StatsCounter nodesMade("Kd-Tree Accelerator", "Interior kd-tree nodes made"); // NOBOOK
		++nodesMade; // NOBOOK
		split = s;
		flags &= ~3;
		flags |= axis;
	}
	float SplitPos() const { return split; }
	int nPrimitives() const { return nPrims >> 2; }
	int SplitAxis() const { return flags & 3; }
	bool IsLeaf() const { return (flags & 3) == 3; }
	union {
		u_int flags;   // Both
		float split;   // Interior
		u_int nPrims;  // Leaf
	};
	union {
		u_int aboveChild;           // Interior
		MailboxPrim *onePrimitive;  // Leaf
		MailboxPrim **primitives;   // Leaf
	};
};
struct BoundEdge {
	// BoundEdge Public Methods
	BoundEdge() { }
	BoundEdge(float tt, int pn, bool starting) {
		t = tt;
		primNum = pn;
		type = starting ? START : END;
	}
	bool operator<(const BoundEdge &e) const {
		if (t == e.t)
			return (int)type < (int)e.type;
		else return t < e.t;
	}
	float t;
	int primNum;
	enum { START, END } type;
};
// KdTreeAccel Declarations
struct KdAccelNode;
class  KdTreeAccel : public Aggregate {
public:
	// KdTreeAccel Public Methods
	KdTreeAccel(const vector<Reference<Primitive> > &p,
		int icost, int scost,
		float ebonus, int maxp, int maxDepth);
	BBox WorldBound() const { return bounds; }
	bool CanIntersect() const { return true; }
	~KdTreeAccel();
	void buildTree(int nodeNum, const BBox &bounds,
	    const vector<BBox> &primBounds,
		int *primNums, int nprims, int depth,
		BoundEdge *edges[3],
		int *prims0, int *prims1, int badRefines = 0);
	bool Intersect(const Ray &ray, Intersection *isect) const;
	bool IntersectP(const Ray &ray) const;
private:
	// KdTreeAccel Private Data
	int isectCost, traversalCost, maxPrims;
	float emptyBonus;
	u_int nMailboxes;
	MailboxPrim *mailboxPrims;
	mutable int curMailboxId;
	KdAccelNode *nodes;
	int nAllocedNodes, nextFreeNode;
	BBox bounds;
	MemoryArena arena;
};
struct KdToDo {
	const KdAccelNode *node;
	float tmin, tmax;
};
// KdTreeAccel Method Definitions
KdTreeAccel::
    KdTreeAccel(const vector<Reference<Primitive> > &p,
		int icost, int tcost,
		float ebonus, int maxp, int maxDepth)
	: isectCost(icost), traversalCost(tcost),
	maxPrims(maxp), emptyBonus(ebonus) {
	vector<Reference<Primitive > > prims;
	for (u_int i = 0; i < p.size(); ++i)
		p[i]->FullyRefine(prims);
	// Initialize mailboxes for _KdTreeAccel_
	curMailboxId = 0;
	nMailboxes = prims.size();
	mailboxPrims = (MailboxPrim *)AllocAligned(nMailboxes *
		sizeof(MailboxPrim));
	for (u_int i = 0; i < nMailboxes; ++i)
		new (&mailboxPrims[i]) MailboxPrim(prims[i]);
	// Build kd-tree for accelerator
	nextFreeNode = nAllocedNodes = 0;
	if (maxDepth <= 0)
		maxDepth =
		    Round2Int(8 + 1.3f * Log2Int(float(prims.size())));
	// Compute bounds for kd-tree construction
	vector<BBox> primBounds;
	primBounds.reserve(prims.size());
	for (u_int i = 0; i < prims.size(); ++i) {
		BBox b = prims[i]->WorldBound();
		bounds = Union(bounds, b);
		primBounds.push_back(b);
	}
	// Allocate working memory for kd-tree construction
	BoundEdge *edges[3];
	for (int i = 0; i < 3; ++i)
		edges[i] = new BoundEdge[2*prims.size()];
	int *prims0 = new int[prims.size()];
	int *prims1 = new int[(maxDepth+1) * prims.size()];
	// Initialize _primNums_ for kd-tree construction
	int *primNums = new int[prims.size()];
	for (u_int i = 0; i < prims.size(); ++i)
		primNums[i] = i;
	// Start recursive construction of kd-tree
	buildTree(0, bounds, primBounds, primNums,
	          prims.size(), maxDepth, edges,
			  prims0, prims1);
	// Free working memory for kd-tree construction
	delete[] primNums;
	for (int i = 0; i < 3; ++i)
		delete[] edges[i];
	delete[] prims0;
	delete[] prims1;
}
KdTreeAccel::~KdTreeAccel() {
	for (u_int i = 0; i < nMailboxes; ++i)
		mailboxPrims[i].~MailboxPrim();
	FreeAligned(mailboxPrims);
	FreeAligned(nodes);
}
void KdTreeAccel::buildTree(int nodeNum,
        const BBox &nodeBounds,
		const vector<BBox> &allPrimBounds, int *primNums,
		int nPrims, int depth, BoundEdge *edges[3],
		int *prims0, int *prims1, int badRefines) {
	Assert(nodeNum == nextFreeNode); // NOBOOK
	// Get next free node from _nodes_ array
	if (nextFreeNode == nAllocedNodes) {
		int nAlloc = max(2 * nAllocedNodes, 512);
		KdAccelNode *n = (KdAccelNode *)AllocAligned(nAlloc *
			sizeof(KdAccelNode));
		if (nAllocedNodes > 0) {
			memcpy(n, nodes,
			       nAllocedNodes * sizeof(KdAccelNode));
			FreeAligned(nodes);
		}
		nodes = n;
		nAllocedNodes = nAlloc;
	}
	++nextFreeNode;
	// Initialize leaf node if termination criteria met
	if (nPrims <= maxPrims || depth == 0) {
		nodes[nodeNum].initLeaf(primNums, nPrims,
		                       mailboxPrims, arena);
		return;
	}
	// Initialize interior node and continue recursion
	// Choose split axis position for interior node
	int bestAxis = -1, bestOffset = -1;
	float bestCost = INFINITY;
	float oldCost = isectCost * float(nPrims);
	Vector d = nodeBounds.pMax - nodeBounds.pMin;
	float totalSA = (2.f * (d.x*d.y + d.x*d.z + d.y*d.z));
	float invTotalSA = 1.f / totalSA;
	// Choose which axis to split along
	int axis;
	if (d.x > d.y && d.x > d.z) axis = 0;
	else axis = (d.y > d.z) ? 1 : 2;
	int retries = 0;
	retrySplit:
	// Initialize edges for _axis_
	for (int i = 0; i < nPrims; ++i) {
		int pn = primNums[i];
		const BBox &bbox = allPrimBounds[pn];
		edges[axis][2*i] =
		    BoundEdge(bbox.pMin[axis], pn, true);
		edges[axis][2*i+1] =
			BoundEdge(bbox.pMax[axis], pn, false);
	}
	sort(&edges[axis][0], &edges[axis][2*nPrims]);
	// Compute cost of all splits for _axis_ to find best
	int nBelow = 0, nAbove = nPrims;
	for (int i = 0; i < 2*nPrims; ++i) {
		if (edges[axis][i].type == BoundEdge::END) --nAbove;
		float edget = edges[axis][i].t;
		if (edget > nodeBounds.pMin[axis] &&
			edget < nodeBounds.pMax[axis]) {
			// Compute cost for split at _i_th edge
			int otherAxis[3][2] = { {1,2}, {0,2}, {0,1} };
			int otherAxis0 = otherAxis[axis][0];
			int otherAxis1 = otherAxis[axis][1];
			float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
			             		(edget - nodeBounds.pMin[axis]) *
				                (d[otherAxis0] + d[otherAxis1]));
			float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
								(nodeBounds.pMax[axis] - edget) *
								(d[otherAxis0] + d[otherAxis1]));
			float pBelow = belowSA * invTotalSA;
			float pAbove = aboveSA * invTotalSA;
			float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.f;
			float cost = traversalCost + isectCost * (1.f - eb) *
				(pBelow * nBelow + pAbove * nAbove);
			// Update best split if this is lowest cost so far
			if (cost < bestCost)  {
				bestCost = cost;
				bestAxis = axis;
				bestOffset = i;
			}
		}
		if (edges[axis][i].type == BoundEdge::START) ++nBelow;
	}
	Assert(nBelow == nPrims && nAbove == 0); // NOBOOK
	// Create leaf if no good splits were found
	if (bestAxis == -1 && retries < 2) {
		++retries;
		axis = (axis+1) % 3;
		goto retrySplit;
	}
	if (bestCost > oldCost) ++badRefines;
	if ((bestCost > 4.f * oldCost && nPrims < 16) ||
		bestAxis == -1 || badRefines == 3) {
		nodes[nodeNum].initLeaf(primNums, nPrims,
		                     mailboxPrims, arena);
		return;
	}
	// Classify primitives with respect to split
	int n0 = 0, n1 = 0;
	for (int i = 0; i < bestOffset; ++i)
		if (edges[bestAxis][i].type == BoundEdge::START)
			prims0[n0++] = edges[bestAxis][i].primNum;
	for (int i = bestOffset+1; i < 2*nPrims; ++i)
		if (edges[bestAxis][i].type == BoundEdge::END)
			prims1[n1++] = edges[bestAxis][i].primNum;
	// Recursively initialize children nodes
	float tsplit = edges[bestAxis][bestOffset].t;
	nodes[nodeNum].initInterior(bestAxis, tsplit);
	BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
	bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tsplit;
	buildTree(nodeNum+1, bounds0,
		allPrimBounds, prims0, n0, depth-1, edges,
		prims0, prims1 + nPrims, badRefines);
	nodes[nodeNum].aboveChild = nextFreeNode;
	buildTree(nodes[nodeNum].aboveChild, bounds1, allPrimBounds,
		prims1, n1, depth-1, edges,
		prims0, prims1 + nPrims, badRefines);
}
bool KdTreeAccel::Intersect(const Ray &ray,
		Intersection *isect) const {
	// Compute initial parametric range of ray inside kd-tree extent
	float tmin, tmax;
	if (!bounds.IntersectP(ray, &tmin, &tmax))
		return false;
	// Prepare to traverse kd-tree for ray
	int rayId = curMailboxId++;
	Vector invDir(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
	#define MAX_TODO 64
	KdToDo todo[MAX_TODO];
	int todoPos = 0;
	// Traverse kd-tree nodes in order for ray
	bool hit = false;
	const KdAccelNode *node = &nodes[0];
	while (node != NULL) {
		// Bail out if we found a hit closer than the current node
		if (ray.maxt < tmin) break;
		//static StatsCounter nodesTraversed("Kd-Tree Accelerator", //NOBOOK
		//	"Number of kd-tree nodes traversed by normal rays"); //NOBOOK
		//++nodesTraversed; //NOBOOK
		if (!node->IsLeaf()) {
			// Process kd-tree interior node
			// Compute parametric distance along ray to split plane
			int axis = node->SplitAxis();
			float tplane = (node->SplitPos() - ray.o[axis]) *
				invDir[axis];
			// Get node children pointers for ray
			const KdAccelNode *firstChild, *secondChild;
			int belowFirst = (ray.o[axis] <  node->SplitPos()) ||
					 (ray.o[axis] == node->SplitPos() && ray.d[axis] >= 0);
			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->aboveChild];
			}
			else {
				firstChild = &nodes[node->aboveChild];
				secondChild = node + 1;
			}
			// Advance to next child node, possibly enqueue other child
			if (tplane > tmax || tplane <= 0)
				node = firstChild;
			else if (tplane < tmin)
				node = secondChild;
			else {
				// Enqueue _secondChild_ in todo list
				todo[todoPos].node = secondChild;
				todo[todoPos].tmin = tplane;
				todo[todoPos].tmax = tmax;
				++todoPos;
				node = firstChild;
				tmax = tplane;
			}
		}
		else {
			// Check for intersections inside leaf node
			u_int nPrimitives = node->nPrimitives();
			if (nPrimitives == 1) {
				MailboxPrim *mp = node->onePrimitive;
				// Check one primitive inside leaf node
				if (mp->lastMailboxId != rayId) {
					mp->lastMailboxId = rayId;
					if (mp->primitive->Intersect(ray, isect))
						hit = true;
				}
			}
			else {
				MailboxPrim **prims = node->primitives;
				for (u_int i = 0; i < nPrimitives; ++i) {
					MailboxPrim *mp = prims[i];
					// Check one primitive inside leaf node
					if (mp->lastMailboxId != rayId) {
						mp->lastMailboxId = rayId;
						if (mp->primitive->Intersect(ray, isect))
							hit = true;
					}
				}
			}
			// Grab next node to process from todo list
			if (todoPos > 0) {
				--todoPos;
				node = todo[todoPos].node;
				tmin = todo[todoPos].tmin;
				tmax = todo[todoPos].tmax;
			}
			else
				break;
		}
	}
	return hit;
}
bool KdTreeAccel::IntersectP(const Ray &ray) const {
	// Compute initial parametric range of ray inside kd-tree extent
	float tmin, tmax;
	if (!bounds.IntersectP(ray, &tmin, &tmax))
		return false;
	// Prepare to traverse kd-tree for ray
	int rayId = curMailboxId++;
	Vector invDir(1.f/ray.d.x, 1.f/ray.d.y, 1.f/ray.d.z);
	#define MAX_TODO 64
	KdToDo todo[MAX_TODO];
	int todoPos = 0;
	const KdAccelNode *node = &nodes[0];
	while (node != NULL) {
		// Update kd-tree shadow ray traversal statistics
		//static StatsCounter nodesTraversed("Kd-Tree Accelerator",
		//	"Number of kd-tree nodes traversed by shadow rays");
		//++nodesTraversed;
		if (node->IsLeaf()) {
			// Check for shadow ray intersections inside leaf node
			u_int nPrimitives = node->nPrimitives();
			if (nPrimitives == 1) {
				MailboxPrim *mp = node->onePrimitive;
				if (mp->lastMailboxId != rayId) {
					mp->lastMailboxId = rayId;
					if (mp->primitive->IntersectP(ray))
						return true;
				}
			}
			else {
				MailboxPrim **prims = node->primitives;
				for (u_int i = 0; i < nPrimitives; ++i) {
					MailboxPrim *mp = prims[i];
					if (mp->lastMailboxId != rayId) {
						mp->lastMailboxId = rayId;
						if (mp->primitive->IntersectP(ray))
							return true;
					}
				}
			}
			// Grab next node to process from todo list
			if (todoPos > 0) {
				--todoPos;
				node = todo[todoPos].node;
				tmin = todo[todoPos].tmin;
				tmax = todo[todoPos].tmax;
			}
			else
				break;
		}
		else {
			// Process kd-tree interior node
			// Compute parametric distance along ray to split plane
			int axis = node->SplitAxis();
			float tplane = (node->SplitPos() - ray.o[axis]) *
				invDir[axis];
			// Get node children pointers for ray
			const KdAccelNode *firstChild, *secondChild;
			int belowFirst = (ray.o[axis] <  node->SplitPos()) ||
					 (ray.o[axis] == node->SplitPos() && ray.d[axis] >= 0);
			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->aboveChild];
			}
			else {
				firstChild = &nodes[node->aboveChild];
				secondChild = node + 1;
			}
			// Advance to next child node, possibly enqueue other child
			if (tplane > tmax || tplane <= 0)
				node = firstChild;
			else if (tplane < tmin)
				node = secondChild;
			else {
				// Enqueue _secondChild_ in todo list
				todo[todoPos].node = secondChild;
				todo[todoPos].tmin = tplane;
				todo[todoPos].tmax = tmax;
				++todoPos;
				node = firstChild;
				tmax = tplane;
			}
		}
	}
	return false;
}
extern "C" DLLEXPORT Primitive *CreateAccelerator(const vector<Reference<Primitive> > &prims,
		const ParamSet &ps) {
	int isectCost = ps.FindOneInt("intersectcost", 80);
	int travCost = ps.FindOneInt("traversalcost", 1);
	float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
	int maxPrims = ps.FindOneInt("maxprims", 1);
	int maxDepth = ps.FindOneInt("maxdepth", -1);
	return new KdTreeAccel(prims, isectCost, travCost,
		emptyBonus, maxPrims, maxDepth);
}
