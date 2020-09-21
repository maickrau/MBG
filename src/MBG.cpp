#include <queue>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdint>
#include <string_view>
#include <phmap.h>
#include "fastqloader.h"
#include "CommonUtils.h"
#include "MBGCommon.h"
#include "MBG.h"
#include "VectorWithDirection.h"
#include "TwobitString.h"
#include "FastHasher.h"

std::string revCompRLE(const std::string& original)
{
	static char mapping[5] { 0, 4, 3, 2, 1 };
	std::string result { original.rbegin(), original.rend() };
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = mapping[(int)result[i]];
	}
	return result;
}

TwobitString revCompRLE(const TwobitString& original)
{
	static char mapping[5] { 0, 4, 3, 2, 1 };
	TwobitString result;
	result.resize(original.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result.set(i, mapping[(int)original.get(original.size()-1-i)]);
	}
	return result;
}

class AdjacentMinimizerList
{
public:
	AdjacentMinimizerList() :
		data(),
		lastHash(0)
	{
	}
	TwobitView getView(size_t coord1, size_t coord2, size_t size) const
	{
		return TwobitView { data[coord1], coord2, coord2+size };
	}
	std::pair<size_t, size_t> addString(std::string_view str, HashType currentHash, HashType previousHash, size_t overlap)
	{
		if (data.size() == 0 || lastHash == 0 || previousHash == 0 || previousHash != lastHash)
		{
			data.emplace_back(str);
			lastHash = currentHash;
			return std::make_pair(data.size()-1, 0);
		}
		assert(overlap < str.size());
		data.back().insert(str.begin() + overlap, str.end());
		lastHash = currentHash;
		assert(data.back().size() >= str.size());
		return std::make_pair(data.size()-1, data.back().size() - str.size());
	}
	AdjacentMinimizerList getReverseComplementStorage() const
	{
		AdjacentMinimizerList result;
		result.data.resize(data.size());
		for (size_t i = 0; i < data.size(); i++)
		{
			result.data[i] = revCompRLE(data[i]);
		}
		return result;
	}
	std::pair<size_t, size_t> getRevCompLocation(size_t coord1, size_t coord2, size_t size) const
	{
		assert(data[coord1].size() >= coord2 + size);
		return std::make_pair(coord1, data[coord1].size() - size - coord2);
	}
private:
	std::vector<TwobitString> data;
	HashType lastHash;
};

class AdjacentLengthList
{
public:
	AdjacentLengthList() :
		sums(),
		counts(),
		lastHash(0)
	{
	}
	std::vector<uint16_t> getData(size_t coord1, size_t coord2, size_t size) const
	{
		std::vector<uint16_t> result;
		result.resize(size, 0);
		assert(sums.size() == counts.size());
		assert(sums[coord1].size() == counts[coord1].size());
		for (size_t i = 0; i < size; i++)
		{
			result[i] = int(double(sums[coord1][coord2+i]) / double(counts[coord1][coord2+i]) + .5);
		}
		return result;
	}
	std::pair<size_t, size_t> addData(const std::vector<uint16_t>& lens, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap)
	{
		assert(end > start);
		assert(end <= lens.size());
		if (sums.size() == 0 || lastHash == 0 || previousHash == 0 || previousHash != lastHash)
		{
			sums.emplace_back(lens.begin() + start, lens.begin() + end);
			counts.emplace_back();
			counts.back().resize(end - start, 1);
			lastHash = currentHash;
			assert(sums.back().size() == counts.back().size());
			return std::make_pair(sums.size()-1, 0);
		}
		assert(overlap < lens.size());
		assert(end > start + overlap);
		sums.back().insert(sums.back().end(), lens.begin() + start + overlap, lens.begin() + end);
		counts.back().resize(counts.back().size() + end - (start + overlap), 1);
		assert(counts.back().size() == sums.back().size());
		lastHash = currentHash;
		assert(sums.back().size() >= end - start);
		return std::make_pair(sums.size()-1, sums.back().size() - (end - start));
	}
	void addCounts(const std::vector<uint16_t>& lens, bool fw, size_t start, size_t end, size_t coord1, size_t coord2)
	{
		assert(coord1 < sums.size());
		assert(sums.size() == counts.size());
		assert(coord2 < sums[coord1].size());
		assert(sums[coord1].size() == counts[coord1].size());
		assert(coord2 + end - start <= sums[coord1].size());
		if (fw)
		{
			for (size_t i = 0; i < end - start; i++)
			{
				if (lens[start+i] < std::numeric_limits<uint16_t>::max() - sums[coord1][coord2+i] && counts[coord1][coord2 + i] < std::numeric_limits<uint8_t>::max() - 1)
				{
					sums[coord1][coord2 + i] += lens[start + i];
					counts[coord1][coord2 + i] += 1;
				}
			}
		}
		else
		{
			for (size_t i = 0; i < end - start; i++)
			{
				if (lens[end - 1 - i] < std::numeric_limits<uint16_t>::max() - sums[coord1][coord2+i] && counts[coord1][coord2 + i] < std::numeric_limits<uint8_t>::max() - 1)
				{
					sums[coord1][coord2 + i] += lens[end - 1 - i];
					counts[coord1][coord2 + i] += 1;
				}
			}
		}
	}
	size_t size() const
	{
		size_t total = 0;
		for (size_t i = 0; i < sums.size(); i++)
		{
			total += sums[i].size();
		}
		return total;
	}
private:
	std::vector<std::vector<uint16_t>> sums;
	std::vector<std::vector<uint8_t>> counts;
	HashType lastHash;
};

class HashList
{
public:
	HashList(size_t kmerSize, bool collapseRunLengths) :
		kmerSize(kmerSize),
		collapseRunLengths(collapseRunLengths)
	{}
	std::vector<size_t> coverage;
	std::vector<uint64_t> fakeFwHashes;
	std::vector<uint64_t> fakeBwHashes;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> sequenceOverlap;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCoverage;
	phmap::flat_hash_map<HashType, std::pair<size_t, bool>> hashToNode;
	size_t numSequenceOverlaps() const
	{
		size_t total = 0;
		for (size_t i = 0; i < sequenceOverlap.size(); i++)
		{
			total += sequenceOverlap[std::make_pair(i, true)].size();
			total += sequenceOverlap[std::make_pair(i, false)].size();
		}
		return total;
	}
	size_t getEdgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		std::tie(from, to) = canon(from, to);
		assert(edgeCoverage.at(from).count(to) == 1);
		return edgeCoverage.at(from).at(to);
	}
	size_t getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		std::tie(from, to) = canon(from, to);
		return sequenceOverlap.at(from).at(to);
	}
	void addSequenceOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap)
	{
		std::tie(from, to) = canon(from, to);
		if (sequenceOverlap[from].count(to) == 1) return;
		sequenceOverlap[from][to] = overlap;
	}
	size_t size() const
	{
		return hashSeqPtr.size();
	}
	std::vector<uint16_t> getHashCharacterLength(size_t index) const
	{
		if (collapseRunLengths)
		{
			std::vector<uint16_t> result;
			result.resize(kmerSize, 1);
			return result;
		}
		return hashCharacterLengths.getData(hashCharacterLengthPtr[index].first, hashCharacterLengthPtr[index].second, kmerSize);
	}
	void addHashCharacterLength(const std::vector<uint16_t>& data, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap)
	{
		if (collapseRunLengths) return;
		hashCharacterLengthPtr.push_back(hashCharacterLengths.addData(data, start, end, currentHash, previousHash, overlap));
	}
	void addHashCharacterLength(const std::vector<uint16_t>& data, bool fw, size_t start, size_t end, size_t node)
	{
		if (collapseRunLengths) return;
		hashCharacterLengths.addCounts(data, fw, start, end, hashCharacterLengthPtr[node].first, hashCharacterLengthPtr[node].second);
	}
	TwobitView getHashSequenceRLE(size_t index) const
	{
		return hashSequences.getView(hashSeqPtr[index].first, hashSeqPtr[index].second, kmerSize);
	}
	TwobitView getRevCompHashSequenceRLE(size_t index) const
	{
		auto pos = hashSequences.getRevCompLocation(hashSeqPtr[index].first, hashSeqPtr[index].second, kmerSize);
		return hashSequencesRevComp.getView(pos.first, pos.second, kmerSize);
	}
	void addHashSequenceRLE(std::string_view seq, HashType currentHash, HashType previousHash, size_t overlap)
	{
		hashSeqPtr.push_back(hashSequences.addString(seq, currentHash, previousHash, overlap));
	}
	void buildReverseCompHashSequences()
	{
		hashSequencesRevComp = hashSequences.getReverseComplementStorage();
	}
	std::pair<size_t, bool> getNodeOrNull(std::string_view sequence) const
	{
		HashType fwHash = hash(sequence);
		auto found = hashToNode.find(fwHash);
		if (found == hashToNode.end())
		{
			return std::pair<size_t, bool> { std::numeric_limits<size_t>::max(), true };
		}
		return found->second;
	}
	std::pair<std::pair<size_t, bool>, HashType> addNode( std::string_view sequence, std::string_view reverse, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap, uint64_t fakeFwHash, uint64_t fakeBwHash)
	{
		HashType fwHash = hash(sequence);
		auto found = hashToNode.find(fwHash);
		if (found != hashToNode.end())
		{
			addHashCharacterLength(sequenceCharacterLength, found->second.second, seqCharLenStart, seqCharLenEnd, found->second.first);
			return std::make_pair(found->second, fwHash);
		}
		assert(found == hashToNode.end());
		HashType bwHash = hash(reverse);
		assert(hashToNode.find(bwHash) == hashToNode.end());
		size_t fwNode = size();
		hashToNode[fwHash] = std::make_pair(fwNode, true);
		hashToNode[bwHash] = std::make_pair(fwNode, false);
		addHashSequenceRLE(sequence, fwHash, previousHash, overlap);
		addHashCharacterLength(sequenceCharacterLength, seqCharLenStart, seqCharLenEnd, fwHash, previousHash, overlap);
		assert(coverage.size() == fwNode);
		coverage.emplace_back(0);
		assert(edgeCoverage.size() == fwNode);
		edgeCoverage.emplace_back();
		assert(sequenceOverlap.size() == fwNode);
		sequenceOverlap.emplace_back();
		assert(fakeFwHashes.size() == fwNode);
		fakeFwHashes.emplace_back(fakeFwHash);
		assert(fakeBwHashes.size() == fwNode);
		fakeBwHashes.emplace_back(fakeBwHash);
		return std::make_pair(std::make_pair(fwNode, true), fwHash);
	}
private:
	AdjacentLengthList hashCharacterLengths;
	std::vector<std::pair<size_t, size_t>> hashCharacterLengthPtr;
	std::vector<size_t> hashCharacterCounts;
	AdjacentMinimizerList hashSequences;
	std::vector<std::pair<size_t, size_t>> hashSeqPtr;
	AdjacentMinimizerList hashSequencesRevComp;
	const size_t kmerSize;
	const bool collapseRunLengths;
};


class LazyString
{
public:
	LazyString(TwobitView first, TwobitView second, size_t overlap) :
	first(first),
	second(second),
	overlap(overlap)
	{
		assert(overlap < first.size());
		assert(overlap < second.size());
	}
	char operator[](size_t index) const
	{
		if (index < first.size()) return first[index];
		assert(index >= first.size());
		size_t secondIndex = index - first.size() + overlap;
		assert(secondIndex < second.size());
		return second[secondIndex];
	}
	std::string_view view(size_t start, size_t kmerSize)
	{
		if (str.size() == 0)
		{
			str.reserve(size());
			str += first.toString();
			str += second.toSubstring(overlap);
		}
		assert(str.size() == size());
		return std::string_view { str.data() + start, kmerSize };
	}
	size_t size() const
	{
		return first.size() + second.size() - overlap;
	}
private:
	TwobitView first;
	TwobitView second;
	const size_t overlap;
	std::string str;
};

class TransitiveCleaner
{
public:
	TransitiveCleaner(size_t kmerSize, const HashList& hashlist)
	{
		transitiveMiddle.resize(hashlist.size());
		getMiddles(kmerSize, hashlist);
	}
	std::vector<std::pair<size_t, bool>> insertMiddles(std::vector<std::pair<size_t, bool>> raw) const
	{
		std::vector<std::pair<size_t, bool>> result;
		while (raw.size() >= 2)
		{
			auto from = raw[raw.size()-2];
			auto to = raw[raw.size()-1];
			if (transitiveMiddle[from].count(to) == 0)
			{
				result.push_back(raw.back());
				raw.pop_back();
				continue;
			}
			raw.pop_back();
			assert(transitiveMiddle[from].at(to).size() > 0);
			raw.insert(raw.end(), transitiveMiddle[from].at(to).begin(), transitiveMiddle[from].at(to).end());
			raw.push_back(to);
		}
		result.push_back(raw.back());
		raw.pop_back();
		std::reverse(result.begin(), result.end());
		return result;
	}
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> newSequenceOverlaps;
private:
	void addMiddles(size_t kmerSize, std::pair<size_t, bool> start, std::pair<size_t, bool> end, LazyString& seq, const HashList& list, const std::unordered_set<size_t>& minimizerPrefixes)
	{
		std::vector<std::pair<size_t, bool>> path;
		std::pair<size_t, bool> old = start;
		size_t oldpos = 0;

		uint64_t fwHash = start.second ? list.fakeFwHashes[start.first] : list.fakeBwHashes[start.first];
		uint64_t bwHash = start.second ? list.fakeBwHashes[start.first] : list.fakeFwHashes[start.first];
		FastHasher fwkmerHasher { kmerSize, fwHash, bwHash };
		assert(minimizerPrefixes.count(fwkmerHasher.hash()) == 1);
		for (size_t i = 1; i < seq.size() - kmerSize; i++)
		{
			fwkmerHasher.addChar(seq[i + kmerSize - 1]);
			fwkmerHasher.removeChar(seq[i - 1]);
			size_t hash = fwkmerHasher.hash();
			if (minimizerPrefixes.count(hash) == 1)
			{
				auto here = list.getNodeOrNull(seq.view(i, kmerSize));
				if (here.first == std::numeric_limits<size_t>::max())
				{
					continue;
				}
				path.push_back(here);
				std::pair<size_t, bool> canonFrom;
				std::pair<size_t, bool> canonTo;
				std::tie(canonFrom, canonTo) = canon(old, here);
				newSequenceOverlaps.emplace_back(canonFrom, canonTo, kmerSize - (i - oldpos));
				old = here;
				oldpos = i;
			}
		}

		if (path.size() > 0)
		{
			assert(old != start);
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(old, end);
			newSequenceOverlaps.emplace_back(canonFrom, canonTo, kmerSize - (seq.size() - kmerSize - oldpos));
			transitiveMiddle[start][end] = path;
		}
		else
		{
			assert(old == start);
		}
	}
	std::unordered_set<size_t> getMinimizerPrefixes(size_t kmerSize, const HashList& hashlist)
	{
		std::unordered_set<size_t> result;
		result.insert(hashlist.fakeFwHashes.begin(), hashlist.fakeFwHashes.end());
		result.insert(hashlist.fakeBwHashes.begin(), hashlist.fakeBwHashes.end());
		return result;
	}
	void getMiddles(size_t kmerSize, const HashList& hashlist)
	{
		std::unordered_set<size_t> minimizerPrefixes = getMinimizerPrefixes(kmerSize, hashlist);
		for (size_t i = 0; i < hashlist.sequenceOverlap.size(); i++)
		{
			std::pair<size_t, bool> fw { i, true };
			if (hashlist.sequenceOverlap[fw].size() > 0)
			{
				TwobitView seq = hashlist.getHashSequenceRLE(i);
				for (auto pair : hashlist.sequenceOverlap[fw])
				{
					assert(pair.first.first >= i);
					TwobitView second { pair.first.second ? hashlist.getHashSequenceRLE(pair.first.first) : hashlist.getRevCompHashSequenceRLE(pair.first.first) };
					LazyString lazy { seq, second, pair.second };
					addMiddles(kmerSize, fw, pair.first, lazy, hashlist, minimizerPrefixes);
				}
			}
			std::pair<size_t, bool> bw { i, false };
			if (hashlist.sequenceOverlap[bw].size() > 0)
			{
				TwobitView seq = hashlist.getRevCompHashSequenceRLE(i);
				for (auto pair : hashlist.sequenceOverlap[bw])
				{
					assert(pair.first.first >= i);
					TwobitView second { pair.first.second ? hashlist.getHashSequenceRLE(pair.first.first) : hashlist.getRevCompHashSequenceRLE(pair.first.first) };
					LazyString lazy { seq, second, pair.second };
					addMiddles(kmerSize, bw, pair.first, lazy, hashlist, minimizerPrefixes);
				}
			}
		}
	}
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, std::vector<std::pair<size_t, bool>>>> transitiveMiddle;
};

std::pair<std::string, std::vector<uint16_t>> runLengthEncode(const std::string& original)
{
	assert(original.size() > 0);
	std::string resultStr;
	std::vector<uint16_t> lens;
	resultStr.reserve(original.size());
	lens.reserve(original.size());
	lens.push_back(1);
	switch(original[0])
	{
		case 'a':
		case 'A':
			resultStr.push_back(1);
			break;
		case 'c':
		case 'C':
			resultStr.push_back(2);
			break;
		case 'g':
		case 'G':
			resultStr.push_back(3);
			break;
		case 't':
		case 'T':
			resultStr.push_back(4);
			break;
		default:
			assert(false);
	}
	for (size_t i = 1; i < original.size(); i++)
	{
		if (original[i] == original[i-1])
		{
			lens.back() += 1;
			continue;
		}
		lens.push_back(1);
		switch(original[i])
		{
			case 'a':
			case 'A':
				resultStr.push_back(1);
				break;
			case 'c':
			case 'C':
				resultStr.push_back(2);
				break;
			case 'g':
			case 'G':
				resultStr.push_back(3);
				break;
			case 't':
			case 'T':
				resultStr.push_back(4);
				break;
			default:
				assert(false);
		}
	}
	return std::make_pair(resultStr, lens);
}

std::pair<std::string, std::vector<uint16_t>> noRunLengthEncode(const std::string& original)
{
	assert(original.size() > 0);
	std::string resultStr;
	resultStr.reserve(original.size());
	for (size_t i = 0; i < original.size(); i++)
	{
		switch(original[i])
		{
			case 'a':
			case 'A':
				resultStr.push_back(1);
				break;
			case 'c':
			case 'C':
				resultStr.push_back(2);
				break;
			case 'g':
			case 'G':
				resultStr.push_back(3);
				break;
			case 't':
			case 'T':
				resultStr.push_back(4);
				break;
			default:
				assert(false);
		}
	}
	std::vector<uint16_t> lens;
	lens.resize(original.size(), 1);
	assert(lens.size() == resultStr.size());
	return std::make_pair(resultStr, lens);
}

std::string strFromRevComp(const std::string& revComp)
{
	std::string result;
	result.reserve(revComp.size());
	for (size_t i = 0; i < revComp.size(); i++)
	{
		result.push_back("-ACGT"[revComp[i]]);
	}
	return result;
}

template <typename F>
void findMinimizerPositions(const std::string& sequence, size_t kmerSize, size_t windowSize, F callback)
{
	if (sequence.size() < kmerSize + windowSize) return;
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	std::deque<std::tuple<size_t, uint64_t, uint64_t, uint64_t>> minimizerOrder;
	minimizerOrder.emplace_back(0, fwkmerHasher.hash(), fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	for (size_t i = 0; i < windowSize-1; i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && std::get<1>(minimizerOrder.back()) > hash) minimizerOrder.pop_back();
		minimizerOrder.emplace_back(i+1, hash, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	}
	auto pos = minimizerOrder.begin();
	while (pos != minimizerOrder.end() && std::get<1>(*pos) == std::get<1>(minimizerOrder.front()))
	{
		callback(std::get<0>(*pos), std::get<2>(*pos), std::get<3>(*pos));
		++pos;
	}
	for (size_t i = windowSize-1; kmerSize+i < sequence.size(); i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		auto oldMinimizer = std::get<1>(minimizerOrder.front());
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && std::get<0>(minimizerOrder.front()) <= i + 1 - windowSize) minimizerOrder.pop_front();
		while (minimizerOrder.size() > 0 && std::get<1>(minimizerOrder.back()) > hash) minimizerOrder.pop_back();
		if (minimizerOrder.size() > 0 && oldMinimizer != std::get<1>(minimizerOrder.front()))
		{
			auto pos = minimizerOrder.begin();
			while (pos != minimizerOrder.end() && std::get<1>(*pos) == std::get<1>(minimizerOrder.front()))
			{
				callback(std::get<0>(*pos), std::get<2>(*pos), std::get<3>(*pos));
				++pos;
			}
		}
		if (minimizerOrder.size() == 0 || hash == std::get<1>(minimizerOrder.front())) callback(i+1, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
		minimizerOrder.emplace_back(i+1, hash, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	}
}

void cleanTransitiveEdges(HashList& result, size_t kmerSize)
{
	TransitiveCleaner cleaner { kmerSize, result };
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> addEdgeCoverage;
	std::vector<std::tuple<std::pair<size_t, bool>, std::pair<size_t, bool>, size_t>> removeEdgeCoverage;

	size_t transitiveEdgesBroken = 0;
	for (size_t node = 0; node < result.edgeCoverage.size(); node++)
	{
		std::pair<size_t, bool> fw { node, true };
		for (auto target : result.edgeCoverage[fw])
		{
			std::vector<std::pair<size_t, bool>> vec;
			vec.push_back(fw);
			vec.push_back(target.first);
			vec = cleaner.insertMiddles(vec);
			if (vec.size() == 2) continue;
			transitiveEdgesBroken += 1;
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(vec[0], vec.back());
			removeEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			for (size_t i = 1; i < vec.size(); i++)
			{
				std::tie(canonFrom, canonTo) = canon(vec[i-1], vec[i]);
				addEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				result.coverage[vec[i].first] += target.second;
			}
		}
		std::pair<size_t, bool> bw { node, false };
		for (auto target : result.edgeCoverage[bw])
		{
			std::vector<std::pair<size_t, bool>> vec;
			vec.push_back(bw);
			vec.push_back(target.first);
			vec = cleaner.insertMiddles(vec);
			if (vec.size() == 2) continue;
			transitiveEdgesBroken += 1;
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(vec[0], vec.back());
			removeEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			for (size_t i = 1; i < vec.size(); i++)
			{
				std::tie(canonFrom, canonTo) = canon(vec[i-1], vec[i]);
				addEdgeCoverage.emplace_back(canonFrom, canonTo, target.second);
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				result.coverage[vec[i].first] += target.second;
			}
		}
	}
	for (auto t : cleaner.newSequenceOverlaps)
	{
		result.sequenceOverlap[std::get<0>(t)][std::get<1>(t)] = std::get<2>(t);
	}
	for (auto t : addEdgeCoverage)
	{
		result.edgeCoverage[std::get<0>(t)][std::get<1>(t)] += std::get<2>(t);
	}
	for (auto t : removeEdgeCoverage)
	{
		assert(result.edgeCoverage.at(std::get<0>(t)).count(std::get<1>(t)) == 1);
		assert(result.edgeCoverage.at(std::get<0>(t)).at(std::get<1>(t)) >= std::get<2>(t));
		result.edgeCoverage[std::get<0>(t)][std::get<1>(t)] -= std::get<2>(t);
	}
	std::cerr << transitiveEdgesBroken << " transitive edges cleaned" << std::endl;
}

HashList loadReadsAsHashes(const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const bool hpc, const bool collapseRunLengths)
{
	HashList result { kmerSize, collapseRunLengths };
	size_t totalNodes = 0;
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&result, &totalNodes, kmerSize, windowSize, hpc](const FastQ& read){
			if (read.sequence.size() == 0) return;
			std::string seq;
			std::vector<uint16_t> lens;
			if (hpc)
			{
				std::tie(seq, lens) = runLengthEncode(read.sequence);
			}
			else
			{
				std::tie(seq, lens) = noRunLengthEncode(read.sequence);
			}
			if (seq.size() <= kmerSize + windowSize) return;
			std::string revSeq = revCompRLE(seq);
			size_t lastMinimizerPosition = std::numeric_limits<size_t>::max();
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			HashType lastHash = 0;
			findMinimizerPositions(seq, kmerSize, windowSize, [kmerSize, windowSize, &lastHash, &last, &lens, &seq, &revSeq, &result, &lastMinimizerPosition, &totalNodes](size_t pos, uint64_t fwHash, uint64_t bwHash)
			{
				assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos > lastMinimizerPosition);
				assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= windowSize);
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				std::string_view minimizerSequence { seq.data() + pos, kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
				std::pair<size_t, bool> current;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				std::tie(current, lastHash) = result.addNode(minimizerSequence, revMinimizerSequence, lens, pos, pos + kmerSize, lastHash, overlap, fwHash, bwHash);
				if (last.first != std::numeric_limits<size_t>::max() && pos - lastMinimizerPosition < kmerSize)
				{
					assert(lastMinimizerPosition + kmerSize >= pos);
					result.addSequenceOverlap(last, current, overlap);
					auto pair = canon(last, current);
					result.edgeCoverage[pair.first][pair.second] += 1;
				}
				lastMinimizerPosition = pos;
				result.coverage[current.first] += 1;
				last = current;
				totalNodes += 1;
			});
		});
	}
	result.buildReverseCompHashSequences();
	std::cerr << totalNodes << " nodes" << std::endl;
	std::cerr << result.size() << " distinct fw/bw sequence nodes" << std::endl;
	return result;
}

class UnitigGraph
{
public:
	std::vector<std::vector<std::pair<NodeType, bool>>> unitigs;
	std::vector<std::vector<size_t>> unitigCoverage;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCov;
	size_t edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw) const
	{
		return edgeCoverage(std::make_pair(from, fromFw), std::make_pair(to, toFw));
	}
	size_t edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		std::tie(from, to) = canon(from, to);
		return edgeCov[from].at(to);
	}
	size_t& edgeCoverage(size_t from, bool fromFw, size_t to, bool toFw)
	{
		return edgeCoverage(std::make_pair(from, fromFw), std::make_pair(to, toFw));
	}
	size_t& edgeCoverage(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
	{
		std::tie(from, to) = canon(from, to);
		return edgeCov[from][to];
	}
	double averageCoverage(size_t i) const
	{
		size_t total = 0;
		for (auto cov : unitigCoverage[i]) total += cov;
		assert(unitigCoverage[i].size() > 0);
		return (double)total / (double)unitigCoverage[i].size();
	}
	UnitigGraph filterNodes(const std::vector<bool>& kept) const
	{
		assert(kept.size() == unitigs.size());
		UnitigGraph result;
		std::vector<size_t> newIndex;
		newIndex.resize(unitigs.size(), std::numeric_limits<size_t>::max());
		size_t nextIndex = 0;
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			if (!kept[i]) continue;
			newIndex[i] = nextIndex;
			nextIndex += 1;
		}
		result.unitigs.resize(nextIndex);
		result.unitigCoverage.resize(nextIndex);
		result.edges.resize(nextIndex);
		result.edgeCov.resize(nextIndex);
		for (size_t i = 0; i < unitigs.size(); i++)
		{
			if (newIndex[i] == std::numeric_limits<size_t>::max()) continue;
			result.unitigs[newIndex[i]] = unitigs[i];
			result.unitigCoverage[newIndex[i]] = unitigCoverage[i];
			std::pair<size_t, bool> fw { i, true };
			std::pair<size_t, bool> bw { i, false };
			std::pair<size_t, bool> newFw { newIndex[i], true };
			std::pair<size_t, bool> newBw { newIndex[i], false };
			for (auto to : edges[fw])
			{
				if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
				result.edges[newFw].emplace(newIndex[to.first], to.second);
			}
			for (auto to : edges[bw])
			{
				if (newIndex[to.first] == std::numeric_limits<size_t>::max()) continue;
				result.edges[newBw].emplace(newIndex[to.first], to.second);
			}
			for (auto pair : edgeCov[fw])
			{
				if (newIndex[pair.first.first] == std::numeric_limits<size_t>::max()) continue;
				result.edgeCov[newFw][std::make_pair(newIndex[pair.first.first], pair.first.second)] = pair.second;
			}
			for (auto pair : edgeCov[bw])
			{
				if (newIndex[pair.first.first] == std::numeric_limits<size_t>::max()) continue;
				result.edgeCov[newBw][std::make_pair(newIndex[pair.first.first], pair.first.second)] = pair.second;
			}
		}
		return result;
	}
	size_t numNodes() const
	{
		return unitigs.size();
	}
	size_t numEdges() const
	{
		size_t result = 0;
		for (size_t i = 0; i < edges.size(); i++)
		{
			std::pair<size_t, bool> fw { i, true };
			std::pair<size_t, bool> bw { i, false };
			for (auto edge : edges[fw])
			{
				auto c = canon(fw, edge);
				if (c.first == fw && c.second == edge) result += 1;
			}
			for (auto edge : edges[bw])
			{
				auto c = canon(bw, edge);
				if (c.first == bw && c.second == edge) result += 1;
			}
		}
		return result;
	}
};

void startUnitig(UnitigGraph& result, const UnitigGraph& old, std::pair<size_t, bool> start, const VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>>& edges, std::vector<std::pair<size_t, bool>>& belongsToUnitig)
{
	size_t currentUnitig = result.unitigs.size();
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	std::pair<size_t, bool> pos = start;
	assert(belongsToUnitig.at(pos.first).first == std::numeric_limits<size_t>::max());
	belongsToUnitig[pos.first] = std::make_pair(currentUnitig, pos.second);
	if (pos.second)
	{
		result.unitigs.back().insert(result.unitigs.back().end(), old.unitigs[pos.first].begin(), old.unitigs[pos.first].end());
		result.unitigCoverage.back().insert(result.unitigCoverage.back().end(), old.unitigCoverage[pos.first].begin(), old.unitigCoverage[pos.first].end());
	}
	else
	{
		for (size_t i = old.unitigs[pos.first].size()-1; i < old.unitigs[pos.first].size(); i--)
		{
			result.unitigs.back().emplace_back(old.unitigs[pos.first][i].first, !old.unitigs[pos.first][i].second);
			result.unitigCoverage.back().emplace_back(old.unitigCoverage[pos.first][i]);
		}
	}
	while (true)
	{
		if (edges.at(pos).size() != 1) break;
		auto newPos = *edges.at(pos).begin();
		auto revPos = std::make_pair(newPos.first, !newPos.second);
		if (edges.at(revPos).size() != 1) break;
		if (newPos == start)
		{
			result.edges[std::make_pair(currentUnitig, true)].emplace(currentUnitig, true);
			result.edgeCoverage(currentUnitig, true, currentUnitig, true) = old.edgeCoverage(pos.first, pos.second, newPos.first, newPos.second);
			break;
		}
		if (belongsToUnitig.at(newPos.first).first != std::numeric_limits<size_t>::max())
		{
			assert(newPos.first == pos.first);
			assert(newPos.second != pos.second);
			assert(belongsToUnitig.at(newPos.first).first == currentUnitig);
			assert(belongsToUnitig.at(newPos.first).second != newPos.second);
			result.edges[std::make_pair(currentUnitig, belongsToUnitig.at(pos.first).second)].emplace(currentUnitig, !belongsToUnitig.at(pos.first).second);
			result.edgeCoverage(currentUnitig, belongsToUnitig.at(pos.first).second, currentUnitig, !belongsToUnitig.at(pos.first).second) = old.edgeCoverage(pos.first, pos.second, newPos.first, newPos.second);
			break;
		}
		pos = newPos;
		assert(belongsToUnitig.at(pos.first).first == std::numeric_limits<size_t>::max());
		belongsToUnitig[pos.first] = std::make_pair(currentUnitig, pos.second);
		if (pos.second)
		{
			result.unitigs.back().insert(result.unitigs.back().end(), old.unitigs[pos.first].begin(), old.unitigs[pos.first].end());
			result.unitigCoverage.back().insert(result.unitigCoverage.back().end(), old.unitigCoverage[pos.first].begin(), old.unitigCoverage[pos.first].end());
		}
		else
		{
			for (size_t i = old.unitigs[pos.first].size()-1; i < old.unitigs[pos.first].size(); i--)
			{
				result.unitigs.back().emplace_back(old.unitigs[pos.first][i].first, !old.unitigs[pos.first][i].second);
				result.unitigCoverage.back().emplace_back(old.unitigCoverage[pos.first][i]);
			}
		}
	}
}

class SparseEdgeContainer
{
public:
	SparseEdgeContainer(size_t size)
	{
		firstEdge.resize(size, std::make_pair(std::numeric_limits<size_t>::max(), false));
	}
	void addEdge(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
	{
		if (firstEdge[from].first == std::numeric_limits<size_t>::max())
		{
			firstEdge[from] = to;
			return;
		}
		if (firstEdge[from] == to) return;
		extraEdges[from].push_back(to);
	}
	std::vector<std::pair<size_t, bool>> operator[](std::pair<size_t, bool> index) const
	{
		return getEdges(index);
	}
	std::vector<std::pair<size_t, bool>> getEdges(std::pair<size_t, bool> from) const
	{
		std::vector<std::pair<size_t, bool>> result;
		if (firstEdge[from].first == std::numeric_limits<size_t>::max()) return result;
		result.push_back(firstEdge[from]);
		auto found = extraEdges.find(from);
		if (found == extraEdges.end()) return result;
		result.insert(result.end(), found->second.begin(), found->second.end());
		return result;
	}
	size_t size() const
	{
		return firstEdge.size();
	}
private:
	VectorWithDirection<std::pair<size_t, bool>> firstEdge;
	phmap::flat_hash_map<std::pair<size_t, bool>, std::vector<std::pair<size_t, bool>>> extraEdges;
};

void startUnitig(UnitigGraph& result, std::pair<size_t, bool> start, const SparseEdgeContainer& edges, std::vector<bool>& belongsToUnitig, const HashList& hashlist, size_t minCoverage)
{
	size_t currentUnitig = result.unitigs.size();
	result.unitigs.emplace_back();
	result.unitigCoverage.emplace_back();
	result.edges.emplace_back();
	result.edgeCov.emplace_back();
	std::pair<size_t, bool> pos = start;
	assert(!belongsToUnitig[pos.first]);
	belongsToUnitig[pos.first] = true;
	result.unitigs.back().emplace_back(pos);
	result.unitigCoverage.back().emplace_back(hashlist.coverage[pos.first]);
	while (true)
	{
		if (edges[pos].size() != 1) break;
		auto newPos = edges[pos][0];
		auto revPos = std::make_pair(newPos.first, !newPos.second);
		if (edges[revPos].size() != 1) break;
		if (newPos == start)
		{
			break;
		}
		if (belongsToUnitig[newPos.first])
		{
			assert(newPos.first == pos.first);
			assert(newPos.second != pos.second);
			break;
		}
		pos = newPos;
		assert(!belongsToUnitig[pos.first]);
		belongsToUnitig[pos.first] = true;
		result.unitigs.back().emplace_back(pos);
		result.unitigCoverage.back().emplace_back(hashlist.coverage[pos.first]);
	}
}

SparseEdgeContainer getCoveredEdges(const HashList& hashlist, size_t minCoverage)
{
	SparseEdgeContainer result { hashlist.coverage.size() };
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto edge : hashlist.edgeCoverage.at(fw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(fw, edge.first);
			result.addEdge(reverse(edge.first), reverse(fw));
		}
		for (auto edge : hashlist.edgeCoverage.at(bw))
		{
			if (edge.second < minCoverage) continue;
			result.addEdge(bw, edge.first);
			result.addEdge(reverse(edge.first), reverse(bw));
		}
	}
	return result;
}

UnitigGraph getUnitigGraph(const HashList& hashlist, size_t minCoverage)
{
	UnitigGraph result;
	std::vector<bool> belongsToUnitig;
	belongsToUnitig.resize(hashlist.coverage.size(), false);
	std::unordered_map<std::pair<size_t, bool>, std::pair<size_t, bool>> unitigTip;
	auto edges = getCoveredEdges(hashlist, minCoverage);
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage[i] < minCoverage) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		auto fwEdges = edges[fw];
		auto bwEdges = edges[bw];
		if (bwEdges.size() != 1)
		{
			if (!belongsToUnitig[i])
			{
				startUnitig(result, fw, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
			for (auto edge : bwEdges)
			{
				if (belongsToUnitig[edge.first]) continue;
				assert(hashlist.coverage[edge.first] >= minCoverage);
				startUnitig(result, edge, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
		}
		if (fwEdges.size() != 1)
		{
			if (!belongsToUnitig[i])
			{
				startUnitig(result, bw, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
			for (auto edge : fwEdges)
			{
				if (belongsToUnitig[edge.first]) continue;
				assert(hashlist.coverage[edge.first] >= minCoverage);
				startUnitig(result, edge, edges, belongsToUnitig, hashlist, minCoverage);
				assert(result.unitigs.size() > 0);
				unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
				unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
			}
		}
	}
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (belongsToUnitig[i]) continue;
		if (hashlist.coverage[i] < minCoverage) continue;
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		auto fwEdges = edges[fw];
		auto bwEdges = edges[bw];
		assert(fwEdges.size() == 1);
		assert(bwEdges.size() == 1);
		startUnitig(result, fw, edges, belongsToUnitig, hashlist, minCoverage);
		assert(result.unitigs.size() > 0);
		assert(result.unitigs.back().back() == reverse(bwEdges[0]));
		unitigTip[result.unitigs.back().back()] = std::make_pair(result.unitigs.size()-1, true);
		unitigTip[reverse(result.unitigs.back()[0])] = std::make_pair(result.unitigs.size()-1, false);
	}
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage[i] < minCoverage) continue;
		assert(belongsToUnitig[i]);
	}
	for (auto tip : unitigTip)
	{
		auto fromNode = tip.first;
		auto fromUnitig = tip.second;
		for (auto edge : edges[fromNode])
		{
			auto toNodeFw = edge;
			auto toNodeRev = reverse(edge);
			assert(unitigTip.count(toNodeRev) == 1);
			auto toUnitig = reverse(unitigTip.at(toNodeRev));
			assert(hashlist.coverage[fromNode.first] >= minCoverage);
			assert(hashlist.coverage[toNodeFw.first] >= minCoverage);
			result.edges[fromUnitig].emplace(toUnitig);
			result.edges[reverse(toUnitig)].emplace(reverse(fromUnitig));
			result.edgeCoverage(fromUnitig, toUnitig) = hashlist.getEdgeCoverage(fromNode, toNodeFw);
		}
	}
	return result;
}

UnitigGraph getNodeGraph(const HashList& hashlist, size_t minCoverage)
{
	UnitigGraph result;
	std::vector<size_t> newIndex;
	newIndex.resize(hashlist.coverage.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < hashlist.coverage.size(); i++)
	{
		if (hashlist.coverage[i] < minCoverage) continue;
		newIndex[i] = result.unitigs.size();
		result.unitigs.emplace_back();
		result.unitigs.back().emplace_back(i, true);
		result.unitigCoverage.emplace_back();
		result.unitigCoverage.back().emplace_back(hashlist.coverage[i]);
	}
	result.edges.resize(result.unitigs.size());
	result.edgeCov.resize(result.unitigs.size());
	for (size_t i = 0; i < hashlist.edgeCoverage.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto pair : hashlist.edgeCoverage[fw])
		{
			if (pair.second < minCoverage) continue;
			assert(hashlist.coverage[i] >= minCoverage);
			assert(hashlist.coverage[pair.first.first] >= minCoverage);
			assert(newIndex[i] != std::numeric_limits<size_t>::max());
			std::pair<size_t, bool> from { newIndex[fw.first], true };
			std::pair<size_t, bool> to { newIndex[pair.first.first], pair.first.second };
			result.edges[from].emplace(to);
			result.edges[reverse(to)].emplace(reverse(from));
			result.edgeCoverage(from, to) = pair.second;
		}
		for (auto pair : hashlist.edgeCoverage[bw])
		{
			if (pair.second < minCoverage) continue;
			assert(hashlist.coverage[i] >= minCoverage);
			assert(hashlist.coverage[pair.first.first] >= minCoverage);
			assert(newIndex[i] != std::numeric_limits<size_t>::max());
			std::pair<size_t, bool> from { newIndex[bw.first], false };
			std::pair<size_t, bool> to { newIndex[pair.first.first], pair.first.second };
			result.edges[from].emplace(to);
			result.edges[reverse(to)].emplace(reverse(from));
			result.edgeCoverage(from, to) = pair.second;
		}
	}
	return result;
}

UnitigGraph getUnitigs(const UnitigGraph& oldgraph)
{
	UnitigGraph result;
	VectorWithDirection<std::unordered_set<std::pair<size_t, bool>>> edges;
	edges.resize(oldgraph.unitigs.size());
	for (size_t i = 0; i < oldgraph.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : oldgraph.edges[fw])
		{
			edges[fw].emplace(to);
			edges[reverse(to)].emplace(reverse(fw));
		}
		for (auto to : oldgraph.edges[bw])
		{
			edges[bw].emplace(to);
			edges[reverse(to)].emplace(reverse(bw));
		}
	}
	std::vector<std::pair<size_t, bool>> belongsToUnitig;
	belongsToUnitig.resize(oldgraph.unitigs.size(), std::make_pair(std::numeric_limits<size_t>::max(), true));
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		auto bw = std::make_pair(node, false);
		if (edges.at(fw).size() != 1)
		{
			for (auto start : edges.at(fw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, bw, edges, belongsToUnitig);
		}
		if (edges.at(bw).size() != 1)
		{
			for (auto start : edges.at(bw))
			{
				if (belongsToUnitig.at(start.first).first != std::numeric_limits<size_t>::max()) continue;
				startUnitig(result, oldgraph, start, edges, belongsToUnitig);
			}
			if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig);
		}
	}
	for (size_t node = 0; node < oldgraph.unitigs.size(); node++)
	{
		auto fw = std::make_pair(node, true);
		if (belongsToUnitig.at(node).first == std::numeric_limits<size_t>::max()) startUnitig(result, oldgraph, fw, edges, belongsToUnitig);
	}
	for (size_t i = 0; i < oldgraph.edges.size(); i++)
	{
		auto fw = std::make_pair(i, true);
		for (auto curr : oldgraph.edges.at(fw))
		{
			auto from = belongsToUnitig.at(fw.first);
			bool prevFw = fw.second;
			auto to = belongsToUnitig.at(curr.first);
			bool currFw = curr.second;
			if (from.first == to.first) continue;
			from.second = !(from.second ^ prevFw);
			to.second = !(to.second ^ currFw);
			result.edges[from].emplace(to);
			result.edgeCoverage(from, to) = oldgraph.edgeCoverage(fw, curr);
		}
		auto bw = std::make_pair(i, false);
		for (auto curr : oldgraph.edges.at(bw))
		{
			auto from = belongsToUnitig.at(bw.first);
			bool prevFw = bw.second;
			auto to = belongsToUnitig.at(curr.first);
			bool currFw = curr.second;
			if (from.first == to.first) continue;
			from.second = !(from.second ^ prevFw);
			to.second = !(to.second ^ currFw);
			result.edges[from].emplace(to);
			result.edgeCoverage(from, to) = oldgraph.edgeCoverage(bw, curr);
		}
	}
	return result;
}

std::string getSequence(const std::string& rle, const std::vector<uint16_t>& characterLength)
{
	std::string result;
	assert(rle.size() == characterLength.size());
	for (size_t i = 0; i < rle.size(); i++)
	{
		for (size_t j = 0; j < characterLength[i]; j++)
		{
			result += "-ACGT"[(int)rle[i]];
		}
	}
	return result;
}

size_t getOverlapFromRLE(const HashList& hashlist, std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	size_t overlap = hashlist.getOverlap(from, to);
	size_t RLEoverlap = 0;
	auto lens = hashlist.getHashCharacterLength(to.first);
	assert(lens.size() > overlap);
	for (size_t offset = 0; offset < overlap; offset++)
	{
		size_t i;
		if (to.second)
		{
			i = offset;
		}
		else
		{
			i = lens.size()-offset-1;
		}
		RLEoverlap += lens[i];
	}
	return RLEoverlap;
}

void writeGraph(const UnitigGraph& unitigs, const std::string& filename, const HashList& hashlist)
{
	std::ofstream file { filename };
	file << "H\tVN:Z:1.0" << std::endl;
	for (size_t i = 0; i < unitigs.unitigs.size(); i++)
	{
		file << "S\t" << i << "\t";
		size_t length = 0;
		for (size_t j = 0; j < unitigs.unitigs[i].size(); j++)
		{
			auto to = unitigs.unitigs[i][j];
			std::string sequenceRLE;
			std::vector<uint16_t> sequenceCharacterLength = hashlist.getHashCharacterLength(to.first);
			if (to.second)
			{
				sequenceRLE = hashlist.getHashSequenceRLE(to.first).toString();
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first).toString();
				std::reverse(sequenceCharacterLength.begin(), sequenceCharacterLength.end());
			}
			if (j > 0)
			{
				auto from = unitigs.unitigs[i][j-1];
				size_t overlap = hashlist.getOverlap(from, to);
				assert(overlap < sequenceRLE.size());
				sequenceRLE = sequenceRLE.substr(overlap);
				sequenceCharacterLength.erase(sequenceCharacterLength.begin(), sequenceCharacterLength.begin() + overlap);
			}
			std::string sequence = getSequence(sequenceRLE, sequenceCharacterLength);
			file << sequence;
			length += sequence.size();
		}
		file << "\tll:f:" << unitigs.averageCoverage(i);
		file << "\tFC:f:" << (unitigs.averageCoverage(i) * length);
		file << std::endl;
	}
	for (size_t i = 0; i < unitigs.edges.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		std::pair<size_t, bool> bw { i, false };
		for (auto to : unitigs.edges[fw])
		{
			std::pair<size_t, bool> last = unitigs.unitigs[i].back();
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t overlap = getOverlapFromRLE(hashlist, last, first);
			file << "L\t" << fw.first << "\t" << (fw.second ? "+" : "-") << "\t" << to.first << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(fw, to) << std::endl;
		}
		for (auto to : unitigs.edges[bw])
		{
			std::pair<size_t, bool> last = reverse(unitigs.unitigs[i][0]);
			std::pair<size_t, bool> first;
			if (to.second)
			{
				first = unitigs.unitigs[to.first][0];
			}
			else
			{
				first = reverse(unitigs.unitigs[to.first].back());
			}
			size_t overlap = getOverlapFromRLE(hashlist, last, first);
			file << "L\t" << bw.first << "\t" << (bw.second ? "+" : "-") << "\t" << to.first << "\t" << (to.second ? "+" : "-") << "\t" << overlap << "M\tec:i:" << unitigs.edgeCoverage(bw, to) << std::endl;
		}
	}
}

auto getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

UnitigGraph filterUnitigsByCoverage(const UnitigGraph& graph, const double filter)
{
	std::vector<bool> kept;
	kept.resize(graph.unitigs.size(), true);
	for (size_t i = 0; i < graph.unitigs.size(); i++)
	{
		if (graph.averageCoverage(i) < filter) kept[i] = false;
	}
	UnitigGraph filtered = graph.filterNodes(kept);
	return filtered;
}

std::pair<size_t, size_t> getSizeAndN50(const HashList& hashlist, const UnitigGraph& graph)
{
	size_t total = 0;
	std::vector<size_t> sizes;
	for (const auto& unitig : graph.unitigs)
	{
		size_t length = 0;
		for (size_t j = 0; j < unitig.size(); j++)
		{
			auto to = unitig[j];
			std::string sequenceRLE;
			std::vector<uint16_t> sequenceCharacterLength = hashlist.getHashCharacterLength(to.first);
			if (to.first)
			{
				sequenceRLE = hashlist.getHashSequenceRLE(to.first).toString();
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first).toString();
				std::reverse(sequenceCharacterLength.begin() ,sequenceCharacterLength.end());
			}
			if (j > 0)
			{
				auto from = unitig[j-1];
				size_t overlap = hashlist.getOverlap(from, to);
				assert(overlap < sequenceRLE.size());
				sequenceRLE = sequenceRLE.substr(overlap);
				sequenceCharacterLength.erase(sequenceCharacterLength.begin(), sequenceCharacterLength.begin() + overlap);
			}
			length += getSequence(sequenceRLE, sequenceCharacterLength).size();
		}
		total += length;
		sizes.push_back(length);
	}
	std::sort(sizes.begin(), sizes.end());
	size_t partialSum = 0;
	for (size_t i = sizes.size()-1; i < sizes.size(); i--)
	{
		partialSum += sizes[i];
		if (partialSum >= total * 0.5) return std::make_pair(total, sizes[i]);
	}
	return std::make_pair(total, 0);
}

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool collapseRunLengths)
{
	auto beforeReading = getTime();
	auto reads = loadReadsAsHashes(inputReads, kmerSize, windowSize, hpc, collapseRunLengths);
	auto beforeCleaning = getTime();
	cleanTransitiveEdges(reads, kmerSize);
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage) unitigs = getUnitigs(filterUnitigsByCoverage(unitigs, minUnitigCoverage));
	auto beforeWrite = getTime();
	writeGraph(unitigs, outputGraph, reads);
	auto beforeStats = getTime();
	auto unitigStats = getSizeAndN50(reads, unitigs);
	auto afterStats = getTime();
	std::cerr << "reading and hashing sequences took " << formatTime(beforeReading, beforeCleaning) << std::endl;
	std::cerr << "cleaning transitive edges took " << formatTime(beforeCleaning, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeWrite) << std::endl;
	std::cerr << "writing the graph took " << formatTime(beforeWrite, beforeStats) << std::endl;
	std::cerr << "calculating stats took " << formatTime(beforeStats, afterStats) << std::endl;
	std::cerr << "nodes: " << unitigs.numNodes() << std::endl;
	std::cerr << "edges: " << unitigs.numEdges() << std::endl;
	std::cerr << "assembly size " << unitigStats.first << " bp, N50 " << unitigStats.second << std::endl;
	std::cerr << "approximate number of k-mers ~ " << (unitigStats.first - unitigs.numNodes() * (kmerSize - windowSize/2 - 1)) << std::endl;
}
