#include <limits>
#include "HPCConsensus.h"
#include "MBGCommon.h"
#include "VectorView.h"
#include "ConsensusMaker.h"

void addCounts(ConsensusMaker& consensusMaker, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const size_t seqStart, const size_t seqEnd, const size_t unitig, const size_t unitigStart, const size_t unitigEnd, const bool fw)
{
	consensusMaker.addStrings(unitig, unitigStart, unitigEnd, [seqStart, seqEnd, &consensusMaker, &seq, &poses, &rawSeq, fw](size_t i)
	{
		size_t seqOff = seqStart + i;
		if (!fw) seqOff = seqEnd - 1 - i;
		assert(seqOff < seq.size());
		uint16_t compressed = seq[seqOff];
		if (!fw) compressed = complement(compressed);
		std::string seq;
		if (fw)
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			seq = rawSeq.substr(expandedStart, expandedEnd - expandedStart);
		}
		else
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			seq = revCompRaw(rawSeq.substr(expandedStart, expandedEnd - expandedStart));
		}
		return std::make_pair(compressed, seq);
	});
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> getHPCUnitigSequences(const HashList& hashlist, const UnitigGraph& unitigs, const std::vector<std::string>& filenames, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	ConsensusMaker consensusMaker;
	std::vector<std::tuple<size_t, size_t, bool>> kmerPosition;
	kmerPosition.resize(hashlist.size(), std::tuple<size_t, size_t, bool> { std::numeric_limits<size_t>::max(), 0, true });
	std::vector<size_t> unitigLengths;
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
		unitigLengths.push_back(offset + kmerSize);
	}
	consensusMaker.init(unitigLengths);
	iterateReadsMultithreaded(filenames, numThreads, [&consensusMaker, &partIterator, &hashlist, &kmerPosition, kmerSize](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&consensusMaker, &hashlist, &kmerPosition, kmerSize](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, uint64_t minHash, const std::vector<size_t>& positions)
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
				VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
				std::pair<size_t, bool> current;
				current = hashlist.getNodeOrNull(minimizerSequence);
				if (current.first == std::numeric_limits<size_t>::max())
				{
					if (currentUnitig != std::numeric_limits<size_t>::max())
					{
						addCounts(consensusMaker, seq, poses, rawSeq, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
					}
					currentUnitig = std::numeric_limits<size_t>::max();
					continue;
				}
				assert(current.first != std::numeric_limits<size_t>::max());
				assert(current.first < kmerPosition.size());
				assert(std::get<0>(kmerPosition[current.first]) != std::numeric_limits<size_t>::max());
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
				addCounts(consensusMaker, seq, poses, rawSeq, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
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
				addCounts(consensusMaker, seq, poses, rawSeq, currentSeqStart, currentSeqEnd, currentUnitig, currentUnitigStart, currentUnitigEnd, currentUnitigForward);
			}
		});
	});
	return consensusMaker.getSequences();
}
