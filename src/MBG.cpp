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

using HashType = unsigned __int128;
using NodeType = size_t;

HashType hash(std::string_view sequence)
{
	size_t half = sequence.size() / 2;
	size_t low = std::hash<std::string_view>{}(std::string_view { sequence.begin(), half });
	size_t high = std::hash<std::string_view>{}(std::string_view { (sequence.begin() + half), sequence.size() - half });
	return (HashType)low + (((HashType)high) << 64);
}

std::ostream& operator<<(std::ostream& os, HashType t)
{
	if (t == 0)
	{
		os << "0";
		return os;
	}
	std::string decimal;
	while (t != 0)
	{
		decimal += "0123456789"[t % 10];
		t /= 10;
	}
	std::reverse(decimal.begin(), decimal.end());
	os << decimal;
	return os;
}

std::istream& operator>>(std::istream& is, HashType& t)
{
	std::string decimal;
	is >> decimal;
	t = 0;
	for (size_t i = 0; i < decimal.size(); i++)
	{
		t *= 10;
		if (decimal[i] >= '0' && decimal[i] <= '9') t += decimal[i]-'0';
	}
	return is;
}

namespace std 
{
	template <> struct hash<HashType>
	{
		size_t operator()(HashType x) const
		{
			return (size_t)x ^ (size_t)(x >> 64);
		}
	};
	template <> struct hash<std::pair<HashType, bool>>
	{
		size_t operator()(std::pair<HashType, bool> x) const
		{
			return hash<HashType>{}(x.first);
		}
	};
	template <> struct hash<std::pair<HashType, HashType>>
	{
		size_t operator()(std::pair<HashType, HashType> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
	template <> struct hash<std::pair<size_t, bool>>
	{
		size_t operator()(std::pair<size_t, bool> x) const
		{
			return (size_t)x.first;
		}
	};
	template <> struct hash<std::pair<size_t, size_t>>
	{
		size_t operator()(std::pair<size_t, size_t> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
}

std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos)
{
	return std::make_pair(pos.first, !pos.second);
}

std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (to.first < from.first)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	if (to.first == from.first && !to.second && !from.second)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	return std::make_pair(from, to);
}

template <typename T>
class VectorWithDirection
{
public:
	void resize(size_t newSize)
	{
		forward.resize(newSize);
		backward.resize(newSize);
	}
	void resize(size_t newSize, T item)
	{
		forward.resize(newSize, item);
		backward.resize(newSize, item);
	}
	size_t size() const
	{
		return forward.size();
	}
	template <typename... Args>
	void emplace_back(Args... args)
	{
		forward.emplace_back(args...);
		backward.emplace_back(args...);
	}
	typename std::vector<T>::reference at(std::pair<size_t, bool> index)
	{
		assert(index.first < forward.size());
		return (*this)[index];
	}
	typename std::vector<T>::const_reference at(std::pair<size_t, bool> index) const
	{
		assert(index.first < forward.size());
		return (*this)[index];
	}
	typename std::vector<T>::reference operator[](std::pair<size_t, bool> index)
	{
		if (index.second) return forward[index.first];
		return backward[index.first];
	}
	typename std::vector<T>::const_reference operator[](std::pair<size_t, bool> index) const
	{
		if (index.second) return forward[index.first];
		return backward[index.first];
	}
	void clear()
	{
		forward.clear();
		backward.clear();
	}
private:
	std::vector<T> forward;
	std::vector<T> backward;
};

class FastHasher
{
public:
	FastHasher(size_t kmerSize) :
	fwHash(0),
	bwHash(0),
	kmerSize(kmerSize % 64)
	{
		precalcRots();
	}
	__attribute__((always_inline))
	void addChar(char c)
	{
		fwHash = rotlone(fwHash) ^ fwAddHash(c);
		bwHash = rotrone(bwHash) ^ bwAddHash(c);
	}
	__attribute__((always_inline))
	void removeChar(char c)
	{
		fwHash ^= fwRemoveHash(c);
		bwHash ^= bwRemoveHash(c);
	}
	__attribute__((always_inline))
	uint64_t hash() const
	{
		return std::min(fwHash, bwHash);
	}
private:
	char complement(char c) const
	{
		static std::vector<char> comp { 0, 4, 3, 2, 1 };
		return comp[c];
	}
	__attribute__((always_inline))
	uint64_t rotlone(uint64_t val) const
	{
		return (val << 1) | (val >> (64-1));
	};
	__attribute__((always_inline))
	uint64_t rotrone(uint64_t val) const
	{
		return (val >> 1) | (val << (64-1));
	};
	__attribute__((always_inline))
	uint64_t rotlk(uint64_t val) const
	{
		return (val << kmerSize) | (val >> (64-kmerSize));
	};
	__attribute__((always_inline))
	uint64_t rotlkmin1(uint64_t val) const
	{
		return (val << (kmerSize-1)) | (val >> (64-(kmerSize-1)));
	};
	void precalcRots()
	{
		for (int i = 0; i < 5; i++)
		{
			fwAdd[i] = charHashes[i];
		}
		for (int i = 0; i < 5; i++)
		{
			fwRemove[i] = rotlk(charHashes[i]);
		}
		for (int i = 0; i < 5; i++)
		{
			bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)]);
		}
		for (int i = 0; i < 5; i++)
		{
			bwRemove[i] = rotrone(charHashes[(int)complement(i)]);
		}
	}
	__attribute__((always_inline))
	uint64_t fwAddHash(char c)
	{
		return fwAdd[(int)c];
	}
	__attribute__((always_inline))
	uint64_t fwRemoveHash(char c)
	{
		return fwRemove[(int)c];
	}
	__attribute__((always_inline))
	uint64_t bwAddHash(char c)
	{
		return bwAdd[(int)c];
	}
	__attribute__((always_inline))
	uint64_t bwRemoveHash(char c)
	{
		return bwRemove[(int)c];
	}
	uint64_t fwAdd[5];
	uint64_t fwRemove[5];
	uint64_t bwAdd[5];
	uint64_t bwRemove[5];
	// https://bioinformatics.stackexchange.com/questions/19/are-there-any-rolling-hash-functions-that-can-hash-a-dna-sequence-and-its-revers
	static constexpr uint64_t hashA = 0x3c8bfbb395c60474;
	static constexpr uint64_t hashC = 0x3193c18562a02b4c;
	static constexpr uint64_t hashG = 0x20323ed082572324;
	static constexpr uint64_t hashT = 0x295549f54be24456;
	uint64_t charHashes[5] { 0, 0x3c8bfbb395c60474, 0x3193c18562a02b4c, 0x20323ed082572324, 0x295549f54be24456 };
	uint64_t fwHash;
	uint64_t bwHash;
	size_t kmerSize;
};

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

class AdjacentMinimizerList
{
public:
	AdjacentMinimizerList() :
		data(),
		lastHash(0)
	{
	}
	std::string_view getView(size_t coord1, size_t coord2, size_t size) const
	{
		return std::string_view { data[coord1].data() + coord2, size };
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
		data.back() += std::string { str.begin() + overlap, str.end() };
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
	std::vector<std::string> data;
	HashType lastHash;
};

class AdjacentLengthList
{
public:
	AdjacentLengthList() :
		data(),
		lastHash(0)
	{
	}
	std::vector<uint16_t> getData(size_t coord1, size_t coord2, size_t size) const
	{
		std::vector<uint16_t> result { data[coord1].begin() + coord2, data[coord1].begin() + coord2 + size };
		return result;
	}
	std::pair<size_t, size_t> addData(const std::vector<uint16_t>& lens, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap)
	{
		assert(end > start);
		assert(end <= lens.size());
		if (data.size() == 0 || lastHash == 0 || previousHash == 0 || previousHash != lastHash)
		{
			data.emplace_back(lens.begin() + start, lens.begin() + end);
			lastHash = currentHash;
			return std::make_pair(data.size()-1, 0);
		}
		assert(overlap < lens.size());
		assert(end > start + overlap);
		data.back().insert(data.back().end(), lens.begin() + start + overlap, lens.begin() + end);
		lastHash = currentHash;
		assert(data.back().size() >= end - start);
		return std::make_pair(data.size()-1, data.back().size() - (end - start));
	}
	size_t size() const
	{
		size_t total = 0;
		for (size_t i = 0; i < data.size(); i++)
		{
			total += data[i].size();
		}
		return total;
	}
private:
	std::vector<std::vector<uint16_t>> data;
	HashType lastHash;
};

class HashList
{
public:
	HashList(size_t kmerSize) :
		kmerSize(kmerSize)
	{}
	std::vector<size_t> coverage;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> sequenceOverlap;
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> edgeCoverage;
	phmap::flat_hash_map<HashType, std::pair<size_t, bool>> hashToNode;
	size_t getOverlap(std::pair<size_t, bool> from, std::pair<size_t, bool> to) const
	{
		std::tie(from, to) = canon(from, to);
		if (sequenceOverlap[from].count(to) == 1) return sequenceOverlap[from].at(to);
		if (sequenceOverlap[from].count(reverse(to)) == 1) return sequenceOverlap[from].at(reverse(to));
		if (sequenceOverlap[reverse(from)].count(to) == 1) return sequenceOverlap[reverse(from)].at(to);
		if (sequenceOverlap[reverse(from)].count(reverse(to)) == 1) return sequenceOverlap[reverse(from)].at(reverse(to));
		return sequenceOverlap[from].at(to);
	}
	size_t size() const
	{
		return hashSeqPtr.size();
	}
	std::vector<uint16_t> getHashCharacterLength(size_t index) const
	{
		return hashCharacterLengths.getData(hashCharacterLengthPtr[index].first, hashCharacterLengthPtr[index].second, kmerSize);
	}
	void addHashCharacterLength(const std::vector<uint16_t>& data, size_t start, size_t end, HashType currentHash, HashType previousHash, size_t overlap)
	{
		hashCharacterLengthPtr.push_back(hashCharacterLengths.addData(data, start, end, currentHash, previousHash, overlap));
	}
	std::string_view getHashSequenceRLE(size_t index) const
	{
		return hashSequences.getView(hashSeqPtr[index].first, hashSeqPtr[index].second, kmerSize);
	}
	std::string_view getRevCompHashSequenceRLE(size_t index) const
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
private:
	AdjacentLengthList hashCharacterLengths;
	std::vector<std::pair<size_t, size_t>> hashCharacterLengthPtr;
	AdjacentMinimizerList hashSequences;
	std::vector<std::pair<size_t, size_t>> hashSeqPtr;
	AdjacentMinimizerList hashSequencesRevComp;
	const size_t kmerSize;
};

std::pair<size_t, bool> getNodeOrNull(const HashList& list, std::string_view sequence)
{
	HashType fwHash = hash(sequence);
	auto found = list.hashToNode.find(fwHash);
	if (found == list.hashToNode.end())
	{
		return std::pair<size_t, bool> { std::numeric_limits<size_t>::max(), true };
	}
	return found->second;
}

class TransitiveCleaner
{
public:
	TransitiveCleaner(size_t kmerSize, const HashList& hashlist)
	{
		transitiveMiddle.resize(hashlist.sequenceOverlap.size());
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
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> newSequenceOverlap;
private:
	void addMiddles(size_t kmerSize, std::pair<size_t, bool> start, std::pair<size_t, bool> end, const std::string& seq, const HashList& list, const std::unordered_set<size_t>& minimizerPrefixes)
	{
		std::vector<std::pair<size_t, bool>> path;
		std::pair<size_t, bool> old = start;
		size_t oldpos = 0;
		size_t kmer = 0;
		for (size_t i = 0; i < 32; i++)
		{
			kmer <<= 2;
			assert(seq[i] >= 1);
			assert(seq[i] <= 5);
			kmer += seq[i]-1;
		}
		for (size_t i = 1; i < seq.size() - kmerSize; i++)
		{
			kmer <<= 2;
			assert(i+32-1 < seq.size());
			assert(seq[i+32-1] >= 1);
			assert(seq[i+32-1] <= 4);
			kmer += seq[i+32-1]-1;
			if (minimizerPrefixes.count(kmer) == 1)
			{
				auto here = getNodeOrNull(list, std::string_view { seq.data() + i, kmerSize });
				if (here.first == std::numeric_limits<size_t>::max()) continue;
				path.push_back(here);
				newSequenceOverlap[old][here] = kmerSize - (i - oldpos);
				newSequenceOverlap[reverse(here)][reverse(old)] = kmerSize - (i - oldpos);
				old = here;
				oldpos = i;
			}
		}
		newSequenceOverlap[old][end] = kmerSize - (seq.size() - kmerSize - oldpos);
		newSequenceOverlap[reverse(end)][reverse(old)] = kmerSize - (seq.size() - kmerSize - oldpos);
		if (path.size() > 0)
		{
			transitiveMiddle[start][end] = path;
		}
	}
	std::unordered_set<size_t> getMinimizerPrefixes(size_t kmerSize, const HashList& hashlist)
	{
		std::unordered_set<size_t> result;
		assert(kmerSize >= 32);
		for (size_t i = 0; i < hashlist.size(); i++)
		{
			size_t kmer = 0;
			for (size_t j = 0; j < 32; j++)
			{
				kmer <<= 2;
				assert(hashlist.getHashSequenceRLE(i)[j] >= 1);
				assert(hashlist.getHashSequenceRLE(i)[j] <= 4);
				kmer += hashlist.getHashSequenceRLE(i)[j]-1;
			}
			result.emplace(kmer);
			kmer = 0;
			for (size_t j = 0; j < 32; j++)
			{
				kmer <<= 2;
				assert(hashlist.getHashSequenceRLE(i)[kmerSize-j-1] >= 1);
				assert(hashlist.getHashSequenceRLE(i)[kmerSize-j-1] <= 4);
				kmer += 3 - (hashlist.getHashSequenceRLE(i)[kmerSize-j-1]-1);
			}
			result.emplace(kmer);
		}
		return result;
	}
	void getMiddles(size_t kmerSize, const HashList& hashlist)
	{
		newSequenceOverlap = hashlist.sequenceOverlap;
		std::unordered_set<size_t> minimizerPrefixes = getMinimizerPrefixes(kmerSize, hashlist);
		std::string seq;
		seq.reserve(kmerSize * 2);
		for (size_t i = 0; i < hashlist.sequenceOverlap.size(); i++)
		{
			std::pair<size_t, bool> fw { i, true };
			if (hashlist.sequenceOverlap[fw].size() > 0)
			{
				seq = hashlist.getHashSequenceRLE(i);
				for (auto pair : hashlist.sequenceOverlap[fw])
				{
					if (pair.first.first < i) continue;
					if (pair.first.second)
					{
						seq += hashlist.getHashSequenceRLE(pair.first.first).substr(pair.second);
					}
					else
					{
						seq += hashlist.getRevCompHashSequenceRLE(pair.first.first).substr(pair.second);
					}
					addMiddles(kmerSize, fw, pair.first, seq, hashlist, minimizerPrefixes);
					seq.resize(kmerSize);
				}
			}
			std::pair<size_t, bool> bw { i, false };
			if (hashlist.sequenceOverlap[bw].size() > 0)
			{
				seq = hashlist.getRevCompHashSequenceRLE(i);
				for (auto pair : hashlist.sequenceOverlap[bw])
				{
					if (pair.first.first < i) continue;
					if (pair.first.second)
					{
						seq += hashlist.getHashSequenceRLE(pair.first.first).substr(pair.second);
					}
					else
					{
						seq += hashlist.getRevCompHashSequenceRLE(pair.first.first).substr(pair.second);
					}
					addMiddles(kmerSize, bw, pair.first, seq, hashlist, minimizerPrefixes);
					seq.resize(kmerSize);
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

void addSequenceOverlap(HashList& list, std::pair<size_t, bool> from, std::pair<size_t, bool> to, const size_t overlap)
{
	if (list.sequenceOverlap[from].count(to) == 1) return;
	std::tie(from, to) = canon(from, to);
	list.sequenceOverlap[from][to] = overlap;
}

std::pair<std::pair<size_t, bool>, HashType> getNode(HashList& list, std::string_view sequence, std::string_view reverse, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap)
{
	HashType fwHash = hash(sequence);
	auto found = list.hashToNode.find(fwHash);
	if (found != list.hashToNode.end())
	{
		return std::make_pair(found->second, fwHash);
	}
	assert(found == list.hashToNode.end());
	HashType bwHash = hash(reverse);
	assert(list.hashToNode.find(bwHash) == list.hashToNode.end());
	size_t fwNode = list.size();
	list.hashToNode[fwHash] = std::make_pair(fwNode, true);
	list.hashToNode[bwHash] = std::make_pair(fwNode, false);
	list.addHashSequenceRLE(sequence, fwHash, previousHash, overlap);
	list.addHashCharacterLength(sequenceCharacterLength, seqCharLenStart, seqCharLenEnd, fwHash, previousHash, overlap);
	assert(list.coverage.size() == fwNode);
	list.coverage.emplace_back(0);
	assert(list.edgeCoverage.size() == fwNode);
	list.edgeCoverage.emplace_back();
	assert(list.sequenceOverlap.size() == fwNode);
	list.sequenceOverlap.emplace_back();
	return std::make_pair(std::make_pair(fwNode, true), fwHash);
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
	std::deque<std::pair<size_t, size_t>> minimizerOrder;
	minimizerOrder.emplace_back(0, fwkmerHasher.hash());
	for (size_t i = 0; i < windowSize-1; i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && minimizerOrder.back().second > hash) minimizerOrder.pop_back();
		minimizerOrder.emplace_back(i+1, hash);
	}
	auto pos = minimizerOrder.begin();
	while (pos != minimizerOrder.end() && pos->second == minimizerOrder.front().second)
	{
		callback(pos->first);
		++pos;
	}
	for (size_t i = windowSize-1; kmerSize+i < sequence.size(); i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		auto oldMinimizer = minimizerOrder.front().second;
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && minimizerOrder.front().first <= i + 1 - windowSize) minimizerOrder.pop_front();
		while (minimizerOrder.size() > 0 && minimizerOrder.back().second > hash) minimizerOrder.pop_back();
		if (minimizerOrder.size() > 0 && oldMinimizer != minimizerOrder.front().second)
		{
			auto pos = minimizerOrder.begin();
			while (pos != minimizerOrder.end() && pos->second == minimizerOrder.front().second)
			{
				callback(pos->first);
				++pos;
			}
		}
		if (minimizerOrder.size() == 0 || hash == minimizerOrder.front().second) callback(i+1);
		minimizerOrder.emplace_back(i+1, hash);
	}
}

void cleanTransitiveEdges(HashList& result, size_t kmerSize)
{
	TransitiveCleaner cleaner { kmerSize, result };
	VectorWithDirection<phmap::flat_hash_map<std::pair<size_t, bool>, size_t>> newEdgeCoverage;

	std::vector<size_t> newCoverage = result.coverage;
	newEdgeCoverage.resize(result.edgeCoverage.size());
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
			if (vec.size() != 2) transitiveEdgesBroken += 1;
			for (size_t i = 1; i < vec.size(); i++)
			{
				newEdgeCoverage[vec[i-1]][vec[i]] += target.second;
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				newCoverage[vec[i].first] += target.second;
			}
		}
		std::pair<size_t, bool> bw { node, false };
		for (auto target : result.edgeCoverage[bw])
		{
			std::vector<std::pair<size_t, bool>> vec;
			vec.push_back(bw);
			vec.push_back(target.first);
			vec = cleaner.insertMiddles(vec);
			if (vec.size() != 2) transitiveEdgesBroken += 1;
			for (size_t i = 1; i < vec.size(); i++)
			{
				newEdgeCoverage[vec[i-1]][vec[i]] += target.second;
			}
			for (size_t i = 1; i < vec.size()-1; i++)
			{
				newCoverage[vec[i].first] += target.second;
			}
		}
	}
	result.sequenceOverlap = cleaner.newSequenceOverlap;
	result.edgeCoverage = newEdgeCoverage;
	result.coverage = newCoverage;
	std::cerr << transitiveEdgesBroken << " transitive edges cleaned" << std::endl;
}

HashList loadReadsAsHashes(const std::string& filename, size_t kmerSize, size_t windowSize, bool hpc)
{
	HashList result { kmerSize };
	size_t totalNodes = 0;
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
		findMinimizerPositions(seq, kmerSize, windowSize, [kmerSize, windowSize, &lastHash, &last, &lens, &seq, &revSeq, &result, &lastMinimizerPosition, &totalNodes](size_t pos)
		{
			assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos > lastMinimizerPosition);
			assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= windowSize);
			assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
			std::string_view minimizerSequence { seq.data() + pos, kmerSize };
			size_t revPos = seq.size() - (pos + kmerSize);
			std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
			std::pair<size_t, bool> current;
			size_t overlap = lastMinimizerPosition + kmerSize - pos;
			std::tie(current, lastHash) = getNode(result, minimizerSequence, revMinimizerSequence, lens, pos, pos + kmerSize, lastHash, overlap);
			if (last.first != std::numeric_limits<size_t>::max() && pos - lastMinimizerPosition < kmerSize)
			{
				assert(lastMinimizerPosition + kmerSize >= pos);
				addSequenceOverlap(result, last, current, overlap);
				auto pair = canon(last, current);
				result.edgeCoverage[pair.first][pair.second] += 1;
			}
			lastMinimizerPosition = pos;
			result.coverage[current.first] += 1;
			last = current;
			totalNodes += 1;
		});
	});
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
				sequenceRLE = hashlist.getHashSequenceRLE(to.first);
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first);
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
				sequenceRLE = hashlist.getHashSequenceRLE(to.first);
			}
			else
			{
				sequenceRLE = hashlist.getRevCompHashSequenceRLE(to.first);
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

int main(int argc, char** argv)
{
	std::string inputReads { argv[1] };
	std::string outputGraph { argv[2] };
	size_t kmerSize = std::stoi(argv[3]);
	size_t windowSize = std::stoi(argv[4]);
	size_t minCoverage = std::stoi(argv[5]);
	double minUnitigCoverage = std::stod(argv[6]);
	bool hpc = true;
	if (argc > 7)
	{
		if (strncmp(argv[7], "1", 1) == 0)
		{
			hpc = false;
		}
	}

	if (windowSize > kmerSize)
	{
		std::cerr << "Window size cannot be greater than k-mer size" << std::endl;
		std::abort();
	}

	auto beforeReading = getTime();
	auto reads = loadReadsAsHashes(inputReads, kmerSize, windowSize, hpc);
	auto beforeCleaning = getTime();
	cleanTransitiveEdges(reads, kmerSize);
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigs(getNodeGraph(reads, minCoverage));
	auto beforeFilter = getTime();
	if (minUnitigCoverage != minCoverage) unitigs = getUnitigs(filterUnitigsByCoverage(unitigs, minUnitigCoverage));
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
	std::cerr << "assembly size " << unitigStats.first << " N50 " << unitigStats.second << std::endl;
}