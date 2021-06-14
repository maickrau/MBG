#include <limits>
#include "HPCConsensus.h"

// allow multiple threads to update the same contig sequence but in different regions
// each mutex covers MutexLength bp in one contig, overlap by half
// ------
//    ------
//       ------
// etc
// size arbitrarily 1Mbp, so ~(3000 + num_contigs) mutexes in a human genome, hopefully not too many
// and the chance of two random hifis falling in the same 1Mbp bucket is ~.3% so hopefully not too much waiting
constexpr size_t MutexLength = 1000000;
constexpr size_t MutexSpan = MutexLength / 2;

void addCounts(std::vector<std::pair<std::string, std::vector<uint16_t>>>& result, std::vector<std::vector<uint8_t>>& runLengthCounts, std::vector<std::vector<bool>>& locked, std::vector<std::mutex*>& seqMutexes, const std::string& seq, const std::vector<uint16_t>& lens, const size_t seqStart, const size_t seqEnd, const size_t unitig, const size_t unitigStart, const size_t unitigEnd, const bool fw)
{
	assert(unitigEnd - unitigStart == seqEnd - seqStart);
	size_t lowMutexIndex = unitigStart / MutexSpan;
	size_t highMutexIndex = (unitigEnd) / MutexSpan;
	std::vector<std::lock_guard<std::mutex>*> guards;
	for (size_t i = lowMutexIndex; i < highMutexIndex; i++)
	{
		guards.emplace_back(new std::lock_guard<std::mutex>{seqMutexes[unitig][i]});
	}
	for (size_t i = 0; i < seqEnd - seqStart; i++)
	{
		size_t off = unitigStart + i;
		if (!fw) off = unitigEnd - 1 - i;
		assert(off < result[unitig].first.size());
		if (locked[unitig][off]) continue;
		if (result[unitig].first[off] == 0)
		{
			if (fw)
			{
				result[unitig].first[off] = seq[seqStart+i];
			}
			else
			{
				result[unitig].first[off] = 5 - seq[seqStart+i];
			}
		}
		else
		{
			assert(!fw || result[unitig].first[off] == seq[seqStart + i]);
			assert(fw || result[unitig].first[off] == 5 - seq[seqStart + i]);
		}
		if (std::numeric_limits<uint16_t>::max() - result[unitig].second[off] < lens[seqStart+i])
		{
			locked[unitig][off] = true;
			continue;
		}
		if (std::numeric_limits<uint8_t>::max() - runLengthCounts[unitig][off] < 1)
		{
			locked[unitig][off] = true;
			continue;
		}
		result[unitig].second[off] += lens[seqStart+i];
		runLengthCounts[unitig][off] += 1;
	}
	for (size_t i = 0; i < guards.size(); i++)
	{
		delete guards[i];
	}
}

std::vector<std::pair<std::string, std::vector<uint16_t>>> getHPCUnitigSequences(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	std::vector<std::pair<std::string, std::vector<uint16_t>>> result;
	std::vector<std::vector<uint8_t>> runLengthCounts;
	std::vector<std::vector<bool>> locked;
	result.resize(unitigs.unitigs.size());
	runLengthCounts.resize(unitigs.unitigs.size());
	locked.resize(unitigs.unitigs.size());
	std::vector<std::tuple<size_t, size_t, bool>> kmerPosition;
	kmerPosition.resize(hashlist.size(), std::tuple<size_t, size_t, bool> { std::numeric_limits<size_t>::max(), 0, true });
	std::vector<std::mutex*> seqMutexes;
	seqMutexes.resize(unitigs.unitigs.size());
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
		locked[i].resize(offset + kmerSize, false);
		rleSize += result[i].first.size();
		size_t numMutexes = 0;
		for (size_t j = 0; j < offset + kmerSize; j += MutexSpan)
		{
			numMutexes += 1;
		}
		assert(numMutexes > 0);
		seqMutexes[i] = new std::mutex[numMutexes];
	}
	iterateReadsMultithreaded(filenames, numThreads, [&result, &locked, &seqMutexes, &runLengthCounts, &partIterator, &hashlist, &kmerPosition, kmerSize](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&result, &locked, &seqMutexes, &runLengthCounts, &hashlist, &kmerPosition, kmerSize](const std::string& seq, const std::vector<uint16_t>& lens, uint64_t minHash, const std::vector<size_t>& positions)
		{
			size_t currentSeqStart = 0;
			size_t currentSeqEnd = 0;
			size_t currentUnitig = std::numeric_limits<size_t>::max();
			size_t currentUnitigStart = 0;
			size_t currentUnitigEnd = 0;
			size_t currentDiagonal = 0;
			bool currentUnitigForward = true;
			std::vector<std::tuple<size_t, size_t, size_t>> matchBlocks;
			for (auto pos : positions)
			{
				std::string_view minimizerSequence { seq.data() + pos, kmerSize };
				std::pair<size_t, bool> current;
				current = hashlist.getNodeOrNull(minimizerSequence);
				assert(current.first != std::numeric_limits<size_t>::max());
				assert(current.first < kmerPosition.size());
				if (std::get<0>(kmerPosition[current.first]) == std::numeric_limits<size_t>::max())
				{
					if (currentUnitig != std::numeric_limits<size_t>::max())
					{
						addCounts(result, runLengthCounts, locked, seqMutexes, seq, lens, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
					}
					currentUnitig = std::numeric_limits<size_t>::max();
				}
				size_t unitig = std::get<0>(kmerPosition[current.first]);
				size_t offset = std::get<1>(kmerPosition[current.first]);
				bool fw = std::get<2>(kmerPosition[current.first]);
				if (!current.second) fw = !fw;
				size_t diagonal;
				if (fw)
				{
					diagonal = pos - offset;
				}
				else
				{
					diagonal = pos + offset;
				}
				if (unitig == currentUnitig && currentUnitigForward == fw && diagonal == currentDiagonal && pos <= currentSeqEnd)
				{
					assert(pos + kmerSize > currentSeqEnd);
					currentSeqEnd = pos + kmerSize;
					if (fw)
					{
						assert(offset + kmerSize > currentUnitigEnd);
						currentUnitigEnd = offset + kmerSize;
					}
					else
					{
						assert(offset < currentUnitigStart);
						currentUnitigStart = offset;
					}
					continue;
				}
				if (currentUnitig == std::numeric_limits<size_t>::max())
				{
					currentUnitig = unitig;
					currentSeqStart = pos;
					currentSeqEnd = pos + kmerSize;
					currentUnitigStart = offset;
					currentUnitigEnd = offset + kmerSize;
					currentDiagonal = diagonal;
					currentUnitigForward = fw;
					continue;
				}
				addCounts(result, runLengthCounts, locked, seqMutexes, seq, lens, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
				currentUnitig = unitig;
				currentSeqStart = pos;
				currentSeqEnd = pos + kmerSize;
				currentUnitigStart = offset;
				currentUnitigEnd = offset + kmerSize;
				currentDiagonal = diagonal;
				currentUnitigForward = fw;
			};
			if (currentUnitig != std::numeric_limits<size_t>::max())
			{
				addCounts(result, runLengthCounts, locked, seqMutexes, seq, lens, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
			}
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
	for (size_t i = 0; i < seqMutexes.size(); i++)
	{
		delete [] seqMutexes[i];
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
