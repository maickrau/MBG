#include "SparseEdgeContainer.h"

SparseEdgeContainer::SparseEdgeContainer(size_t size)
{
	firstEdge.resize(size, std::numeric_limits<uint32_t>::max());
}

void SparseEdgeContainer::addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (pairToInt(to) != std::numeric_limits<uint32_t>::max())
	{
		if (firstEdge[from] == std::numeric_limits<uint32_t>::max())
		{
			firstEdge[from] = pairToInt(to);
			return;
		}
		if (firstEdge[from] == pairToInt(to)) return;
	}
	extraEdges[from].push_back(to);
}

std::vector<std::pair<size_t, bool>> SparseEdgeContainer::operator[](std::pair<size_t, bool> index) const
{
	return getEdges(index);
}

std::vector<std::pair<size_t, bool>> SparseEdgeContainer::getEdges(std::pair<size_t, bool> from) const
{
	std::vector<std::pair<size_t, bool>> result;
	if (firstEdge[from] != std::numeric_limits<uint32_t>::max())
	{
		result.push_back(intToPair(firstEdge[from]));
	}
	auto found = extraEdges.find(from);
	if (found == extraEdges.end()) return result;
	result.insert(result.end(), found->second.begin(), found->second.end());
	return result;
}

size_t SparseEdgeContainer::size() const
{
	return firstEdge.size();
}

uint32_t SparseEdgeContainer::pairToInt(std::pair<size_t, bool> value) const
{
	if (value.first >= (size_t)std::numeric_limits<uint32_t>::max() / 2) return std::numeric_limits<uint32_t>::max();
	return (uint32_t)value.first * 2 + (value.second ? 1 : 0);
}

std::pair<size_t, bool> SparseEdgeContainer::intToPair(uint32_t value) const
{
	return std::make_pair(value / 2, (value % 2) == 1);
}
