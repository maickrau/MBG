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
				size_t lowMutexIndex = offset / MutexSpan;
				size_t highMutexIndex = (offset + kmerSize) / MutexSpan;
				std::vector<std::lock_guard<std::mutex>*> guards;
				for (size_t i = lowMutexIndex; i < highMutexIndex; i++)
				{
					guards.emplace_back(new std::lock_guard<std::mutex>{seqMutexes[unitig][i]});
				}
				for (size_t i = 0; i < kmerSize; i++)
				{
					size_t off = offset + i;
					if (!fw) off = offset + kmerSize - 1 - i;
					if (locked[unitig][off]) continue;
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
					if (std::numeric_limits<uint16_t>::max() - result[unitig].second[off] < lens[pos+i])
					{
						locked[unitig][off] = true;
						continue;
					}
					if (std::numeric_limits<uint8_t>::max() - runLengthCounts[unitig][off] < 1)
					{
						locked[unitig][off] = true;
						continue;
					}
					result[unitig].second[off] += lens[pos+i];
					runLengthCounts[unitig][off] += 1;
				}
				for (size_t i = 0; i < guards.size(); i++)
				{
					delete guards[i];
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
