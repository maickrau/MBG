#ifndef SparseEdgeContainer_h
#define SparseEdgeContainer_h

#include <tuple>
#include <vector>
#include <phmap.h>
#include "VectorWithDirection.h"

class SparseEdgeContainer
{
public:
	SparseEdgeContainer(size_t size);
	void addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	std::vector<std::pair<size_t, bool>> operator[](std::pair<size_t, bool> index) const;
	std::vector<std::pair<size_t, bool>> getEdges(std::pair<size_t, bool> from) const;
	size_t size() const;
private:
	VectorWithDirection<std::pair<size_t, bool>> firstEdge;
	phmap::flat_hash_map<std::pair<size_t, bool>, std::vector<std::pair<size_t, bool>>> extraEdges;
};

#endif
