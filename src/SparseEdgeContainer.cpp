#include "SparseEdgeContainer.h"

SparseEdgeContainer::SparseEdgeContainer(size_t size)
{
	firstEdge.resize(size, std::make_pair(std::numeric_limits<size_t>::max(), false));
}

void SparseEdgeContainer::addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (firstEdge[from].first == std::numeric_limits<size_t>::max())
	{
		firstEdge[from] = to;
		return;
	}
	if (firstEdge[from] == to) return;
	extraEdges[from].push_back(to);
}

std::vector<std::pair<size_t, bool>> SparseEdgeContainer::operator[](std::pair<size_t, bool> index) const
{
	return getEdges(index);
}

std::vector<std::pair<size_t, bool>> SparseEdgeContainer::getEdges(std::pair<size_t, bool> from) const
{
	std::vector<std::pair<size_t, bool>> result;
	if (firstEdge[from].first == std::numeric_limits<size_t>::max()) return result;
	result.push_back(firstEdge[from]);
	auto found = extraEdges.find(from);
	if (found == extraEdges.end()) return result;
	result.insert(result.end(), found->second.begin(), found->second.end());
	return result;
}

size_t SparseEdgeContainer::size() const
{
	return firstEdge.size();
}
