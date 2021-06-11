#include <limits>
#include "HPCConsensus.h"

std::vector<std::pair<std::string, std::vector<uint16_t>>> getHPCUnitigSequences(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const size_t kmerSize, const ReadpartIterator& partIterator)
{
	std::vector<std::pair<std::string, std::vector<uint16_t>>> result;
	std::vector<std::vector<uint8_t>> runLengthCounts;
	result.resize(unitigs.unitigs.size());
	runLengthCounts.resize(unitigs.unitigs.size());
	std::vector<std::tuple<size_t, size_t, bool>> kmerPosition;
	kmerPosition.resize(hashlist.size(), std::tuple<size_t, size_t, bool> { std::numeric_limits<size_t>::max(), 0, true });
	size_t rleSize = 0;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t offset = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			if (j > 0)
			{
				size_t rleOverlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
				assert(rleOverlap < kmerSize);
				offset += kmerSize - rleOverlap;
			}
			assert(std::get<0>(kmerPosition[unitigs.unitigs[i][j].first]) == std::numeric_limits<size_t>::max());
			kmerPosition[unitigs.unitigs[i][j].first] = std::make_tuple(i, offset, unitigs.unitigs[i][j].second);
		}
		result[i].first.resize(offset + kmerSize, 0);
		result[i].second.resize(offset + kmerSize, 0);
		runLengthCounts[i].resize(offset + kmerSize, 0);
		rleSize += result[i].first.size();
	}
	iterateReadsMultithreaded(filenames, 1, [&result, &runLengthCounts, &partIterator, &hashlist, &kmerPosition, kmerSize](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &runLengthCounts, &hashlist, &kmerPosition, kmerSize](const std::string& seq, const std::vector<uint16_t>& lens, uint64_t minHash, const std::vector<size_t>& positions)
		{
			std::string revSeq = revCompRLE(seq);
			for (auto pos : positions)
			{
				std::string_view minimizerSequence { seq.data() + pos, kmerSize };
				std::pair<size_t, bool> current;
				current = hashlist.getNodeOrNull(minimizerSequence);
				assert(current.first != std::numeric_limits<size_t>::max());
				assert(current.first < kmerPosition.size());
				if (std::get<0>(kmerPosition[current.first]) == std::numeric_limits<size_t>::max()) continue;
				size_t unitig = std::get<0>(kmerPosition[current.first]);
				size_t offset = std::get<1>(kmerPosition[current.first]);
				bool fw = std::get<2>(kmerPosition[current.first]);
				if (!current.second) fw = !fw;
				for (size_t i = 0; i < kmerSize; i++)
				{
					size_t off = offset + i;
					if (!fw) off = offset + kmerSize - 1 - i;
					if (result[unitig].first[off] == 0)
					{
						if (fw)
						{
							result[unitig].first[off] = seq[pos+i];
						}
						else
						{
							result[unitig].first[off] = 5 - seq[pos+i];
						}
					}
					else
					{
						assert(!fw || result[unitig].first[off] == seq[pos + i]);
						assert(fw || result[unitig].first[off] == 5 - seq[pos + i]);
					}
					result[unitig].second[off] += lens[pos+i];
					runLengthCounts[unitig][off] += 1;
				}
			};
		});
	});
	size_t expandedSize = 0;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		for (size_t j = 0; j < result[i].second.size(); j++)
		{
			assert(runLengthCounts[i][j] > 0);
			result[i].second[j] = (result[i].second[j] + runLengthCounts[i][j] / 2) / runLengthCounts[i][j];
			expandedSize += result[i].second[j];
		}
	}
	return result;
}

std::string getHPCExpanded(const std::string& seq, const std::vector<uint16_t>& length)
{
	assert(seq.size() == length.size());
	std::string result;
	for (size_t i = 0; i < seq.size(); i++)
	{
		assert(length[i] > 0);
		for (size_t j = 0; j < length[i]; j++)
		{
			switch(seq[i])
			{
				case 1:
					result += 'A';
					break;
				case 2:
					result += 'C';
					break;
				case 3:
					result += 'G';
					break;
				case 4:
					result += 'T';
					break;
				default:
					assert(false);
					break;
			}
		}
	}
	return result;
}
