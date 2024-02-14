#include <limits>
#include <unordered_map>
#include "HPCConsensus.h"
#include "MBGCommon.h"
#include "VectorView.h"
#include "ConsensusMaker.h"
#include "KmerMatcher.h"

using namespace MBG;

namespace MBG
{

void addCounts(ConsensusMaker& consensusMaker, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const size_t seqStart, const size_t seqEnd, const size_t unitig, const size_t unitigStart, const size_t unitigEnd, const bool fw)
{
	consensusMaker.addStrings(unitig, unitigStart, unitigEnd, [seqStart, seqEnd, &consensusMaker, &seq, &poses, &rawSeq, fw](size_t i)
	{
		size_t seqOff = seqStart + i;
		if (!fw) seqOff = seqEnd - 1 - i;
		assert(seqOff < seq.size());
		uint16_t compressed = seq[seqOff];
		if (!fw) compressed = complement(compressed);
		std::variant<size_t, std::string> expanded;
		if (fw)
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			if (compressed >= 0 && compressed <= 3)
			{
				expanded = expandedEnd - expandedStart;
			}
			else
			{
				expanded = rawSeq.substr(expandedStart, expandedEnd - expandedStart);
			}
		}
		else
		{
			size_t expandedStart = poses[seqOff];
			size_t expandedEnd = poses[seqOff+1];
			assert(expandedEnd > expandedStart);
			if (compressed >= 0 && compressed <= 3)
			{
				expanded = expandedEnd - expandedStart;
			}
			else
			{
				expanded = revCompRaw(rawSeq.substr(expandedStart, expandedEnd - expandedStart));
			}
		}
		return std::make_pair(compressed, expanded);
	});
}

std::vector<std::tuple<size_t, size_t, size_t, size_t, bool, size_t, size_t>> getMatchBlocks(const std::vector<std::vector<size_t>>& bpOffsets, const UnitigGraph& unitigs, const ReadPath& path, const std::vector<size_t> unitigLengths, const size_t kmerSize, const size_t readi)
{
	std::vector<std::tuple<size_t, size_t, size_t, size_t, bool, size_t, size_t>> result;
	size_t pathKmerCount = 0;
	std::vector<size_t> pathKmerStarts;
	std::vector<size_t> pathKmerEnds;
	for (size_t i = 0; i < path.path.size(); i++)
	{
		if (i > 0) pathKmerCount -= unitigs.edgeOverlap(path.path[i-1], path.path[i]);
		pathKmerStarts.push_back(pathKmerCount);
		pathKmerCount += unitigs.unitigs[path.path[i].id()].size();
		pathKmerEnds.push_back(pathKmerCount);
	}
	assert(pathKmerStarts[0] == 0);
	assert(pathKmerEnds.back() == pathKmerCount);
	assert(pathKmerCount == path.readPoses.size() + path.leftClip + path.rightClip);
	for (size_t i = 0; i < path.path.size(); i++)
	{
		size_t wantedStart = pathKmerStarts[i];
		size_t wantedEnd = pathKmerEnds[i];
		size_t possibleStart = path.leftClip;
		size_t possibleEnd = pathKmerCount - path.rightClip;
		size_t startKmerIndex = 0;
		size_t endKmerIndex = 0;
		size_t startPosIndex = 0;
		size_t endPosIndex = 0;
		if (wantedStart < possibleStart)
		{
			startPosIndex = 0;
			startKmerIndex = possibleStart - wantedStart;
		}
		else
		{
			startPosIndex = wantedStart - possibleStart;
			startKmerIndex = 0;
		}
		assert(wantedEnd > wantedStart + startKmerIndex);
		assert(possibleEnd > possibleStart + startPosIndex);
		size_t maxMatch = std::min(wantedEnd - wantedStart - startKmerIndex, possibleEnd - possibleStart - startPosIndex);
		endKmerIndex = startKmerIndex + maxMatch - 1;
		endPosIndex = startPosIndex + maxMatch - 1;
		if (!path.path[i].forward())
		{
			startKmerIndex = bpOffsets[path.path[i].id()].size() - 1 - startKmerIndex;
			endKmerIndex = bpOffsets[path.path[i].id()].size() - 1 - endKmerIndex;
			std::swap(startKmerIndex, endKmerIndex);
		}
		assert(i != 0 || startPosIndex == 0);
		assert(i != path.path.size()-1 || endPosIndex == path.readPoses.size()-1);
		assert(startKmerIndex <= endKmerIndex);
		assert(startPosIndex <= endPosIndex);
		assert(endPosIndex < path.readPoses.size());
		assert(endKmerIndex < bpOffsets[path.path[i].id()].size());
		size_t readStart = path.readPoses[startPosIndex];
		size_t readEnd = path.readPoses[endPosIndex] + kmerSize;
		assert(readEnd <= path.readLengthHPC);
		size_t unitigStart = bpOffsets[path.path[i].id()][startKmerIndex];
		size_t unitigEnd = bpOffsets[path.path[i].id()][endKmerIndex] + kmerSize;
		size_t unitig = path.path[i].id();
		bool fw = path.path[i].forward();
		if (unitigStart < unitigs.leftClip[unitig])
		{
			if (fw)
			{
				readStart += unitigs.leftClip[unitig] - unitigStart;
			}
			else
			{
				readEnd -= unitigs.leftClip[unitig] - unitigStart;
			}
			unitigStart = 0;
		}
		else
		{
			unitigStart -= unitigs.leftClip[unitig];
		}
		assert(readEnd <= path.readLengthHPC);
		assert(unitigEnd >= unitigs.leftClip[unitig]);
		unitigEnd -= unitigs.leftClip[unitig];
		if (unitigEnd > unitigLengths[unitig])
		{
			if (fw)
			{
				readEnd -= unitigEnd - unitigLengths[unitig];
			}
			else
			{
				readStart += unitigEnd - unitigLengths[unitig];
			}
			unitigEnd = unitigLengths[unitig];
		}
		assert(readEnd <= path.readLengthHPC);
		assert(readEnd > readStart);
		assert(unitigEnd > unitigStart);
		assert(readEnd - readStart == unitigEnd - unitigStart);
		assert(readEnd <= path.readLengthHPC);
		result.emplace_back(unitig, readStart, unitigStart, readEnd - readStart, fw, i == 0 ? readi : std::numeric_limits<size_t>::max(), i == path.path.size()-1 ? readi : std::numeric_limits<size_t>::max());
	}
	return result;
}

void initializeHelpers(ConsensusMaker& consensusMaker, std::vector<size_t>& unitigLengths, std::vector<std::vector<size_t>>& bpOffsets, const HashList& hashlist, const UnitigGraph& unitigs, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
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
		assert(offset + kmerSize > unitigs.leftClip[i] + unitigs.rightClip[i]);
		unitigLengths.push_back(offset + kmerSize - unitigs.leftClip[i] - unitigs.rightClip[i]);
	}
	consensusMaker.init(unitigLengths);
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::pair<size_t, bool> pos { i, true };
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[edge.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.prepareEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
		pos = std::make_pair(i, false);
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[pos.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.prepareEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
	}
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		std::pair<size_t, bool> pos { i, true };
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[edge.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.addEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
		pos = std::make_pair(i, false);
		for (auto edge : unitigs.edges[pos])
		{
			size_t kmerOverlap = unitigs.edgeOverlap(pos, edge);
			size_t bpOverlap = 0;
			if (kmerOverlap == 0)
			{
				std::pair<size_t, bool> fromKmer = unitigs.unitigs[pos.first].back();
				if (!pos.second) fromKmer = reverse(unitigs.unitigs[pos.first][0]);
				std::pair<size_t, bool> toKmer = unitigs.unitigs[edge.first][0];
				if (!edge.second) toKmer = reverse(unitigs.unitigs[edge.first].back());
				bpOverlap = hashlist.getOverlap(fromKmer, toKmer);
			}
			else
			{
				std::pair<size_t, bool> lastKmer;
				for (size_t i = 0; i < kmerOverlap; i++)
				{
					std::pair<size_t, bool> kmer = unitigs.unitigs[edge.first][i];
					if (!edge.second) kmer = reverse(unitigs.unitigs[edge.first][unitigs.unitigs[edge.first].size() - 1 - i]);
					bpOverlap += kmerSize;
					if (i > 0) bpOverlap -= hashlist.getOverlap(lastKmer, kmer);
					lastKmer = kmer;
				}
			}
			size_t fromClip = pos.second ? unitigs.rightClip[pos.first] : unitigs.leftClip[pos.first];
			size_t toClip = edge.second ? unitigs.leftClip[edge.first] : unitigs.rightClip[edge.first];
			assert(bpOverlap > fromClip + toClip);
			consensusMaker.addEdgeOverlap(pos, edge, bpOverlap - fromClip - toClip);
		}
	}
	consensusMaker.findParentLinks();
}

std::pair<std::vector<CompressedSequenceType>, StringIndex> getHPCUnitigSequences(const HashList& hashlist, const UnitigGraph& unitigs, std::vector<ReadPath>& readPaths, const size_t kmerSize, const ReadpartIterator& partIterator, const size_t numThreads)
{
	ConsensusMaker consensusMaker;
	std::vector<size_t> unitigLengths;
	std::vector<std::vector<size_t>> bpOffsets;
	initializeHelpers(consensusMaker, unitigLengths, bpOffsets, hashlist, unitigs, kmerSize, partIterator, numThreads);
	std::mutex expandedPosMutex;
	std::unordered_map<ReadName, std::vector<size_t>> pathsPerRead;
	for (size_t i = 0; i < readPaths.size(); i++)
	{
		pathsPerRead[readPaths[i].readName].push_back(i);
	}
	partIterator.iterateParts([&consensusMaker, &readPaths, &expandedPosMutex, &hashlist, &unitigLengths, &unitigs, &bpOffsets, &pathsPerRead, kmerSize](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq)
	{
		if (pathsPerRead.count(read.readName) == 0) return;
		std::vector<std::tuple<size_t, size_t, size_t, size_t, bool, size_t, size_t>> matchBlocks;
		for (size_t pathi : pathsPerRead.at(read.readName))
		{
			auto add = getMatchBlocks(bpOffsets, unitigs, readPaths[pathi], unitigLengths, kmerSize, pathi);
			matchBlocks.insert(matchBlocks.end(), add.begin(), add.end());
		}
		for (auto block : matchBlocks)
		{
			size_t unitig = std::get<0>(block);
			size_t readStartPos = std::get<1>(block);
			size_t unitigStartPos = std::get<2>(block);
			size_t matchLength = std::get<3>(block);
			bool forward = std::get<4>(block);
			addCounts(consensusMaker, seq, poses, rawSeq, readStartPos, readStartPos + matchLength, unitig, unitigStartPos, unitigStartPos + matchLength, forward);
			if (std::get<5>(block) != std::numeric_limits<size_t>::max() || std::get<6>(block) != std::numeric_limits<size_t>::max())
			{
				assert(readStartPos + matchLength < poses.size());
				std::lock_guard<std::mutex> lock { expandedPosMutex };
				if (std::get<5>(block) != std::numeric_limits<size_t>::max())
				{
					readPaths[std::get<5>(block)].expandedReadPosStart = poses[readStartPos];
				}
				if (std::get<6>(block) != std::numeric_limits<size_t>::max())
				{
					readPaths[std::get<6>(block)].expandedReadPosEnd = poses[readStartPos + matchLength];
				}
			}
		}
	});
	return consensusMaker.getSequences();
}

void getHpcVariantsAndReadPaths(const HashList& hashlist, const UnitigGraph& unitigs, const size_t kmerSize, ReadpartIterator& partIterator, const size_t numThreads, const double minUnitigCoverage, const size_t minVariantCoverage)
{
	ConsensusMaker consensusMaker;
	std::vector<size_t> unitigLengths;
	std::vector<std::vector<size_t>> bpOffsets;
	initializeHelpers(consensusMaker, unitigLengths, bpOffsets, hashlist, unitigs, kmerSize, partIterator, numThreads);
	std::vector<bool> checkUnitig;
	checkUnitig.resize(unitigs.unitigs.size(), false);
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		double sum = 0;
		for (size_t j = 0; j < unitigs.unitigCoverage[i].size(); j++)
		{
			sum += unitigs.unitigCoverage[i][j];
		}
		double count = unitigs.unitigCoverage[i].size();
		double coverage = sum / count;
		if (coverage >= minUnitigCoverage) checkUnitig[i] = true;
	}
	std::vector<std::tuple<size_t, size_t, bool>> kmerLocator = getKmerLocator(unitigs);
	partIterator.iterateHashes([&consensusMaker, &checkUnitig, &bpOffsets, &unitigs, &hashlist, &unitigLengths, &kmerLocator, kmerSize](const ReadInfo& read, const SequenceCharType& seq, const SequenceLengthType& poses, const std::string& rawSeq, const std::vector<size_t>& positions, const std::vector<HashType>& hashes)
	{
		iterateReadPaths(unitigs, hashlist, kmerSize, kmerLocator, read, positions, hashes, [&consensusMaker, &bpOffsets, &unitigs, &unitigLengths, &checkUnitig, &rawSeq, &poses, &seq, kmerSize](ReadPath path)
		{
			auto blocks = getMatchBlocks(bpOffsets, unitigs, path, unitigLengths, kmerSize, 0);
			for (auto& block : blocks)
			{
				size_t unitig = std::get<0>(block);
				if (!checkUnitig[unitig]) continue;
				size_t readStartPos = std::get<1>(block);
				size_t unitigStartPos = std::get<2>(block);
				size_t matchLength = std::get<3>(block);
				bool forward = std::get<4>(block);
				addCounts(consensusMaker, seq, poses, rawSeq, readStartPos, readStartPos + matchLength, unitig, unitigStartPos, unitigStartPos + matchLength, forward);
			}
		});
	});
	std::vector<HashType> nodeToHash;
	nodeToHash.resize(hashlist.hashToNode.size(), 0);
	for (auto pair : hashlist.hashToNode)
	{
		assert(nodeToHash[pair.second] == 0);
		assert(pair.first != 0);
		nodeToHash[pair.second] = pair.first;
		assert(nodeToHash[pair.second] != 0);
	}
	consensusMaker.prepareHpcVariants(checkUnitig);
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		if (!checkUnitig[i]) continue;
		auto hpcVariants = consensusMaker.getHpcVariants(i, minVariantCoverage);
		for (size_t j = 0; j < hpcVariants.size(); j++)
		{
			size_t variantPos = hpcVariants[j].first;
			size_t kmerIndex = 0;
			while (bpOffsets[i][kmerIndex] + kmerSize <= variantPos)
			{
				kmerIndex += 1;
				assert(kmerIndex < bpOffsets[i].size());
			}
			assert(kmerIndex < bpOffsets[i].size());
			assert(bpOffsets[i][kmerIndex] <= variantPos);
			while (kmerIndex < bpOffsets[i].size() && bpOffsets[i][kmerIndex] <= variantPos)
			{
				bool fw = unitigs.unitigs[i][kmerIndex].second;
				size_t offset = variantPos - bpOffsets[i][kmerIndex];
				assert(offset < kmerSize);
				if (!fw) offset = kmerSize - 1 - offset;
				assert(offset < kmerSize);
				partIterator.addHpcVariants(nodeToHash[unitigs.unitigs[i][kmerIndex].first], offset, hpcVariants[j].second);
				kmerIndex += 1;
			}
		}
	}
}

}
