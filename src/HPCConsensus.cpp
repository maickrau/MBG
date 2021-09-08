#include <limits>
#include "HPCConsensus.h"
#include "MBGCommon.h"
#include "VectorView.h"
#include "ConsensusMaker.h"

class KmerLocator
{
public:
	KmerLocator(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs)
	{
		nodes.emplace_back();
		nodes.back().parent = 0;
		nodes.back().depth = 0;
		std::get<0>(nodes.back().uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<1>(nodes.back().uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<2>(nodes.back().uniqueLocation) = true;
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			for (size_t j = 0; j < unitigs[i].size(); j++)
			{
				addPath(unitigs, i, j, true);
				addPath(unitigs, i, j, false);
			}
		}
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			for (size_t j = 0; j < unitigs[i].size(); j++)
			{
				auto found = locate(unitigs[i], j);
				if (std::get<0>(found) == std::numeric_limits<size_t>::max()) continue;
				assert(std::get<0>(found) == i);
				assert(std::get<1>(found) == j);
				assert(std::get<2>(found) == true);
			}
			auto rev = revCompPath(unitigs[i]);
			for (size_t j = 0; j < rev.size(); j++)
			{
				auto found = locate(rev, j);
				if (std::get<0>(found) == std::numeric_limits<size_t>::max()) continue;
				assert(std::get<0>(found) == i);
				assert(std::get<1>(found) == unitigs[i].size()-1-j);
				assert(std::get<2>(found) == false);
			}
		}
	}
	std::tuple<size_t, size_t, bool> locate(const std::vector<std::pair<size_t, bool>>& path, size_t start) const
	{
		std::pair<size_t, bool> pos = path[start];
		size_t node = 0;
		size_t depth = 0;
		assert(nodes[0].children.count(pos) == 1);
		while (true)
		{
			if (nodes[node].children.count(pos) == 0)
			{
				return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			}
			node = nodes[node].children.at(pos);
			if (std::get<0>(nodes[node].uniqueLocation) != std::numeric_limits<size_t>::max())
			{
				size_t i = std::get<0>(nodes[node].uniqueLocation);
				size_t j = std::get<1>(nodes[node].uniqueLocation);
				bool fw = std::get<2>(nodes[node].uniqueLocation);
				return std::make_tuple(i, j, fw);
			}
			depth += 1;
			if (start+depth >= path.size()) return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
			pos = path[start+depth];
		}
		assert(false);
		return std::make_tuple(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max(), true);
	}
private:
	void addPath(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs, size_t i, size_t j, bool fw)
	{
		size_t node = 0;
		size_t depth = 0;
		while (true)
		{
			size_t checkIndex = j;
			if (fw)
			{
				checkIndex += depth;
			}
			else
			{
				checkIndex -= depth;
			}
			if (checkIndex >= unitigs[i].size()) return;
			std::pair<size_t, bool> checkKmer = unitigs[i][checkIndex];
			if (!fw) checkKmer = reverse(checkKmer);
			if (nodes[node].children.count(checkKmer) == 0) break;
			node = nodes[node].children.at(checkKmer);
			depth += 1;
			if (std::get<0>(nodes[node].uniqueLocation) != std::numeric_limits<size_t>::max()) knockdown(unitigs, node);
		}
		size_t checkIndex = j;
		if (fw)
		{
			checkIndex += depth;
		}
		else
		{
			checkIndex -= depth;
		}
		assert(checkIndex < unitigs[i].size());
		size_t newNode = nodes.size();
		nodes.emplace_back();
		nodes[newNode].parent = node;
		nodes[newNode].depth = depth;
		nodes[newNode].uniqueLocation = std::make_tuple(i, j, fw);
		std::pair<size_t, bool> checkKmer = unitigs[i][checkIndex];
		if (!fw) checkKmer = reverse(checkKmer);
		nodes[node].children[checkKmer] = newNode;
	}
	void knockdown(const std::vector<std::vector<std::pair<size_t, bool>>>& unitigs, size_t node)
	{
		size_t i = std::get<0>(nodes[node].uniqueLocation);
		size_t j = std::get<1>(nodes[node].uniqueLocation);
		bool fw = std::get<2>(nodes[node].uniqueLocation);
		assert(i != std::numeric_limits<size_t>::max());
		size_t depth = nodes[node].depth;
		size_t next;
		if (fw)
		{
			next = j + depth + 1;
		}
		else
		{
			next = j - depth - 1;
		}
		std::get<0>(nodes[node].uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<1>(nodes[node].uniqueLocation) = std::numeric_limits<size_t>::max();
		std::get<2>(nodes[node].uniqueLocation) = true;
		assert(j < unitigs[i].size());
		if (next >= unitigs[i].size()) return;
		std::pair<size_t, bool> key;
		key = unitigs[i][next];
		if (!fw) key = reverse(key);
		assert(nodes[node].children.count(key) == 0);
		size_t nextIndex = nodes.size();
		nodes[node].children[key] = nextIndex;
		nodes.emplace_back();
		nodes[nextIndex].parent = node;
		nodes[nextIndex].depth = depth+1;
		std::get<0>(nodes[nextIndex].uniqueLocation) = i;
		std::get<1>(nodes[nextIndex].uniqueLocation) = j;
		std::get<2>(nodes[nextIndex].uniqueLocation) = fw;
	}
	class TreeNode
	{
	public:
		size_t parent;
		size_t depth;
		std::unordered_map<std::pair<size_t, bool>, size_t> children;
		std::tuple<size_t, size_t, bool> uniqueLocation;
	};
	std::vector<TreeNode> nodes;
};

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
	KmerLocator locator { unitigs.unitigs };
	std::vector<size_t> unitigLengths;
	std::vector<std::vector<size_t>> bpOffsets;
	bpOffsets.resize(unitigs.unitigs.size());
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		size_t offset = 0;
		bpOffsets[i].push_back(0);
		for (size_t j = 1; j < unitigs.unitigs[i].size(); j++)
		{
			size_t rleOverlap = hashlist.getOverlap(unitigs.unitigs[i][j-1], unitigs.unitigs[i][j]);
			assert(rleOverlap < kmerSize);
			offset += kmerSize - rleOverlap;
			bpOffsets[i].push_back(offset);
		}
		unitigLengths.push_back(offset + kmerSize);
	}
	consensusMaker.init(unitigLengths);
	iterateReadsMultithreaded(filenames, numThreads, [&consensusMaker, &unitigLengths, &bpOffsets, &unitigs, &partIterator, &hashlist, &locator, kmerSize](size_t thread, FastQ& read)
	{
		partIterator.iteratePartKmers(read, [&consensusMaker, &hashlist, &unitigLengths, &unitigs, &bpOffsets, &locator, kmerSize](const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, uint64_t minHash, const std::vector<size_t>& positions)
		{
			std::vector<std::pair<size_t, bool>> kmerPath;
			std::vector<size_t> kmerPositions;
			for (auto pos : positions)
			{
				VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
				std::pair<size_t, bool> current = hashlist.getNodeOrNull(minimizerSequence);
				if (current.first == std::numeric_limits<size_t>::max()) continue;
				kmerPath.push_back(current);
				kmerPositions.push_back(pos);
			}
			size_t lastMatchStart = 0;
			size_t lastMatchEnd = 0;
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, bool>> matches;
			for (size_t i = 0; i < kmerPath.size(); i++)
			{
				auto kmerpos = locator.locate(kmerPath, i);
				if (std::get<0>(kmerpos) == std::numeric_limits<size_t>::max()) continue;
				size_t unitig = std::get<0>(kmerpos);
				size_t kmeroffset = std::get<1>(kmerpos);
				bool fw = std::get<2>(kmerpos);
				size_t matchLength = 0;
				while (i + matchLength < kmerPath.size())
				{
					std::pair<size_t, bool> kmer = kmerPath[i + matchLength];
					size_t comparePos;
					if (fw)
					{
						comparePos = kmeroffset + matchLength;
					}
					else
					{
						comparePos = kmeroffset - matchLength;
					}
					if (comparePos >= unitigs.unitigs[unitig].size()) break;
					std::pair<size_t, bool> compareKmer = unitigs.unitigs[unitig][comparePos];
					if (!fw) compareKmer = reverse(compareKmer);
					if (compareKmer != kmer) break;
					if (matchLength > 0)
					{
						size_t seqPosOffset = kmerPositions[i + matchLength] - kmerPositions[i + matchLength - 1];
						size_t unitigPosOffset;
						if (fw)
						{
							assert(comparePos > 0);
							unitigPosOffset = bpOffsets[unitig][comparePos] - bpOffsets[unitig][comparePos - 1];
						}
						else
						{
							assert(comparePos+1 < bpOffsets[unitig].size());
							unitigPosOffset = bpOffsets[unitig][comparePos + 1] - bpOffsets[unitig][comparePos];
						}
						if (seqPosOffset != unitigPosOffset) break;
					}
					matchLength += 1;
				}
				size_t matchPreLength = 0;
				while (i - matchPreLength < kmerPath.size())
				{
					std::pair<size_t, bool> kmer = kmerPath[i - matchPreLength];
					size_t comparePos;
					if (fw)
					{
						comparePos = kmeroffset - matchPreLength;
					}
					else
					{
						comparePos = kmeroffset + matchPreLength;
					}
					if (comparePos >= unitigs.unitigs[unitig].size()) break;
					std::pair<size_t, bool> compareKmer = unitigs.unitigs[unitig][comparePos];
					if (!fw) compareKmer = reverse(compareKmer);
					if (compareKmer != kmer) break;
					if (matchPreLength > 0)
					{
						size_t seqPosOffset = kmerPositions[i - matchPreLength + 1] - kmerPositions[i - matchPreLength];
						size_t unitigPosOffset;
						if (fw)
						{
							assert(comparePos+1 < bpOffsets[unitig].size());
							unitigPosOffset = bpOffsets[unitig][comparePos + 1] - bpOffsets[unitig][comparePos];
						}
						else
						{
							assert(comparePos > 0);
							unitigPosOffset = bpOffsets[unitig][comparePos] - bpOffsets[unitig][comparePos - 1];
						}
						if (seqPosOffset != unitigPosOffset) break;
					}
					matchPreLength += 1;
				}
				assert(matchPreLength >= 1);
				matchPreLength -= 1;
				assert(i - matchPreLength >= lastMatchStart);
				assert(i + matchLength >= lastMatchEnd);
				if (i + matchLength == lastMatchEnd) continue;
				assert(i + matchLength > lastMatchEnd);
				assert(i + matchLength <= kmerPositions.size());
				matches.emplace_back(unitig, i, matchPreLength, matchLength, kmeroffset, fw);
			}
			if (matches.size() == 0) return;
			std::sort(matches.begin(), matches.end(), [](std::tuple<size_t, size_t, size_t, size_t, size_t, bool> left, std::tuple<size_t, size_t, size_t, size_t, size_t, bool> right) { return std::get<1>(left) + std::get<3>(left) < std::get<1>(right) + std::get<3>(right); });
			std::vector<std::tuple<size_t, size_t, size_t, size_t, size_t, bool>> nonOverlappingMatches;
			nonOverlappingMatches.push_back(matches[0]);
			for (size_t i = 1; i < matches.size(); i++)
			{
				assert(nonOverlappingMatches.size() >= 1);
				if (std::get<1>(matches[i]) + std::get<3>(matches[i]) > std::get<1>(nonOverlappingMatches.back()) + std::get<3>(nonOverlappingMatches.back()))
				{
					nonOverlappingMatches.push_back(matches[i]);
					continue;
				}
				assert(std::get<1>(matches[i]) + std::get<3>(matches[i]) == std::get<1>(nonOverlappingMatches.back()) + std::get<3>(nonOverlappingMatches.back()));
				if (std::get<1>(matches[i]) - std::get<2>(matches[i]) < std::get<1>(nonOverlappingMatches.back()) - std::get<2>(nonOverlappingMatches.back()))
				{
					nonOverlappingMatches.pop_back();
					nonOverlappingMatches.push_back(matches[i]);
					continue;
				}
				if (std::get<1>(matches[i]) - std::get<2>(matches[i]) > std::get<1>(nonOverlappingMatches.back()) - std::get<2>(nonOverlappingMatches.back())) continue;
				assert(std::get<1>(matches[i]) - std::get<2>(matches[i]) == std::get<1>(nonOverlappingMatches.back()) - std::get<2>(nonOverlappingMatches.back()));
				assert(std::get<0>(matches[i]) == std::get<0>(nonOverlappingMatches.back()));
			}
			for (size_t i = 0; i < nonOverlappingMatches.size(); i++)
			{
				size_t unitig = std::get<0>(nonOverlappingMatches[i]);
				size_t offset = std::get<1>(nonOverlappingMatches[i]);
				size_t matchPreLength = std::get<2>(nonOverlappingMatches[i]);
				size_t matchLength = std::get<3>(nonOverlappingMatches[i]);
				size_t kmeroffset = std::get<4>(nonOverlappingMatches[i]);
				bool fw = std::get<5>(nonOverlappingMatches[i]);
				assert(matchLength >= 1);
				size_t unitigStart;
				size_t unitigEnd;
				assert(kmeroffset < bpOffsets[unitig].size());
				size_t kmerStart;
				size_t kmerEnd;
				if (fw)
				{
					assert(kmeroffset + matchLength - 1 < bpOffsets[unitig].size());
					assert(matchPreLength <= kmeroffset);
					unitigStart = bpOffsets[unitig][kmeroffset - matchPreLength];
					unitigEnd = bpOffsets[unitig][kmeroffset + matchLength - 1] + kmerSize;
					kmerStart = kmeroffset - matchPreLength;
					kmerEnd = kmeroffset + matchLength;
				}
				else
				{
					assert(kmeroffset - matchLength + 1 < bpOffsets[unitig].size());
					assert(kmeroffset + matchPreLength < bpOffsets[unitig].size());
					unitigStart = bpOffsets[unitig][kmeroffset - matchLength + 1];
					unitigEnd = bpOffsets[unitig][kmeroffset + matchPreLength] + kmerSize;
					kmerStart = kmeroffset - matchLength + 1;
					kmerEnd = kmeroffset + matchPreLength + 1;
				}
				for (size_t j = kmerStart; j < kmerEnd; j++)
				{
					// unitigs.unitigCoverage[unitig][j] += 1;
				}
				assert(unitigStart < unitigEnd);
				assert(unitigEnd <= unitigLengths[unitig]);
				assert(offset + matchLength <= kmerPositions.size());
				assert(matchPreLength <= offset);
				size_t seqStart = kmerPositions[offset - matchPreLength];
				size_t seqEnd = kmerPositions[offset + matchLength - 1] + kmerSize;
				assert(seqEnd > seqStart);
				assert(seqEnd <= seq.size());
				assert(unitigEnd - unitigStart == seqEnd - seqStart);
				addCounts(consensusMaker, seq, poses, rawSeq, seqStart, seqEnd, unitig, unitigStart, unitigEnd, fw);
			};
		});
	});
	return consensusMaker.getSequences();
}
