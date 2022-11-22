#ifndef SparseEdgeContainer_h
#define SparseEdgeContainer_h

#include <tuple>
#include <vector>
#include <phmap.h>
#include "VectorWithDirection.h"

class SparseEdgeContainer
{
public:
	SparseEdgeContainer();
	SparseEdgeContainer(size_t size);
	void addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
	std::vector<std::pair<size_t, bool>> operator[](std::pair<size_t, bool> index) const;
	std::vector<std::pair<size_t, bool>> getEdges(std::pair<size_t, bool> from) const;
	bool hasEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const;
	size_t size() const;
	void emplace_back();
	void resize(size_t newSize);
private:
	uint32_t pairToInt(std::pair<size_t, bool> value) const;
	std::pair<size_t, bool> intToPair(uint32_t value) const;
	VectorWithDirection<uint32_t> firstEdge;
	phmap::flat_hash_map<std::pair<size_t, bool>, std::vector<std::pair<size_t, bool>>> extraEdges;
};

#endif
