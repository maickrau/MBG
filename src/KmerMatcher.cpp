#include "KmerMatcher.h"

using namespace MBG;

namespace MBG
{

std::vector<std::tuple<size_t, size_t, bool>> getKmerLocator(const UnitigGraph& graph)
{
	size_t maxKmer = 0;
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			maxKmer = std::max(maxKmer, graph.unitigs[i][j].first);
		}
	}
	std::vector<std::tuple<size_t, size_t, bool>> kmerLocator;
	kmerLocator.resize(maxKmer+1, std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true));
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		for (size_t j = 0; j < graph.unitigs[i].size(); j++)
		{
			assert(std::get<0>(kmerLocator[graph.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			if (graph.unitigs[i][j].second)
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, j, true);
			}
			else
			{
				kmerLocator[graph.unitigs[i][j].first] = std::make_tuple(i, graph.unitigs[i].size()-1-j, false);
			}
		}
	}
	return kmerLocator;
}

}
