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
#include "MBG.h"

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
	FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash) :
	fwHash(fwHash),
	bwHash(bwHash),
	kmerSize(kmerSize % 64)
	{
		precalcRots();
	}
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
	__attribute__((always_inline))
	uint64_t getFwHash() const
	{
		return fwHash;
	}
	__attribute__((always_inline))
	uint64_t getBwHash() const
	{
		return bwHash;
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

class TwobitString
{
public:
	TwobitString() :
	data(),
	realSize(0)
	{}
	TwobitString(const std::string_view& str) :
	data(),
	realSize(0)
	{
		*this = str;
	}
	char get(size_t i) const
	{
		size_t pos = i / 4;
		size_t off = (i % 4) * 2;
		return ((data[pos] >> off) & 3) + 1;
	}
	void set(size_t i, char c)
	{
		assert(c >= 1);
		assert(c <= 4);
		size_t pos = i / 4;
		size_t off = (i % 4) * 2;
		data[pos] &= ~(3 << off);
		data[pos] |= (c-1) << off;
	}
	size_t size() const
	{
		return realSize;
	}
	TwobitString& operator=(const std::string& str)
	{
		data.clear();
		resize(str.size());
		for (size_t i = 0; i < str.size(); i++)
		{
			set(i, str[i]);
		}
		return *this;
	}
	TwobitString& operator=(const std::string_view& str)
	{
		data.clear();
		resize(str.size());
		for (size_t i = 0; i < str.size(); i++)
		{
			set(i, str[i]);
		}
		return *this;
	}
	void resize(size_t size)
	{
		size_t max = (size + 3) / 4;
		data.resize(max, 0);
		realSize = size;
	}
	void push_back(char c)
	{
		size_t pos = realSize / 4;
		if (pos >= data.size())
		{
			data.resize(data.size() * 2, 0);
		}
		set(realSize, c);
		realSize += 1;
	}
	template <typename Iter>
	void insert(Iter start, Iter end)
	{
		Iter i = start;
		while (i != end)
		{
			push_back(*i);
			++i;
		}
	}
private:
	std::vector<unsigned char> data;
	size_t realSize;
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

class TwobitView
{
public:
	TwobitView(const TwobitString& str, const size_t start, const size_t end) :
	str(str),
	start(start),
	end(end)
	{}
	char operator[](size_t i) const
	{
		return str.get(start+i);
	}
	size_t size() const
	{
		return end-start;
	}
	std::string toString() const
	{
		std::string result;
		result.reserve(end - start);
		for (size_t i = start; i < end; i++)
		{
			result.push_back(str.get(i));
		}
		return result;
	}
	std::string toSubstring(size_t substrStart) const
	{
		assert(start + substrStart < end);
		std::string result;
		result.reserve(end - (start + substrStart));
		for (size_t i = start + substrStart; i < end; i++)
		{
			result.push_back(str.get(i));
		}
		return result;
	}
private:
	const TwobitString& str;
	const size_t start;
	const size_t end;
};

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

std::pair<std::pair<size_t, bool>, HashType> addNode(HashList& list, std::string_view hashableSequence, std::string_view realSequence, std::string_view reverseHashableSequence, std::string_view reverseRealSequence, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap, uint64_t fakeFwHash, uint64_t fakeBwHash)
{
	HashType fwHash = hash(hashableSequence);
	auto found = list.hashToNode.find(fwHash);
	if (found != list.hashToNode.end())
	{
		list.addHashCharacterLength(sequenceCharacterLength, found->second.second, seqCharLenStart, seqCharLenEnd, found->second.first);
		return std::make_pair(found->second, fwHash);
	}
	assert(found == list.hashToNode.end());
	HashType bwHash = hash(reverseHashableSequence);
	assert(list.hashToNode.find(bwHash) == list.hashToNode.end());
	size_t fwNode = list.size();
	list.hashToNode[fwHash] = std::make_pair(fwNode, true);
	list.hashToNode[bwHash] = std::make_pair(fwNode, false);
	list.addHashSequenceRLE(realSequence, fwHash, previousHash, overlap);
	list.addHashCharacterLength(sequenceCharacterLength, seqCharLenStart, seqCharLenEnd, fwHash, previousHash, overlap);
	assert(list.coverage.size() == fwNode);
	list.coverage.emplace_back(0);
	assert(list.edgeCoverage.size() == fwNode);
	list.edgeCoverage.emplace_back();
	assert(list.sequenceOverlap.size() == fwNode);
	list.sequenceOverlap.emplace_back();
	assert(list.fakeFwHashes.size() == fwNode);
	list.fakeFwHashes.emplace_back(fakeFwHash);
	assert(list.fakeBwHashes.size() == fwNode);
	list.fakeBwHashes.emplace_back(fakeBwHash);
	return std::make_pair(std::make_pair(fwNode, true), fwHash);
}

std::pair<std::pair<size_t, bool>, HashType> addNode(HashList& list, std::string_view hash1, std::string_view hash2, std::string_view hash3, std::string_view hash4, std::string_view realSequence, std::string_view revhash1, std::string_view revhash2, std::string_view revhash3, std::string_view revhash4, std::string_view reverseRealSequence, const std::vector<uint16_t>& sequenceCharacterLength, size_t seqCharLenStart, size_t seqCharLenEnd, HashType previousHash, size_t overlap, uint64_t fakeFwHash, uint64_t fakeBwHash)
{
	HashType fwHash1 = hash(hash1);
	HashType fwHash2 = hash(hash2);
	HashType fwHash3 = hash(hash3);
	HashType fwHash4 = hash(hash4);
	auto found1 = list.hashToNode.find(fwHash1);
	auto found2 = list.hashToNode.find(fwHash2);
	auto found3 = list.hashToNode.find(fwHash3);
	auto found4 = list.hashToNode.find(fwHash4);
	HashType bwHash1 = hash(revhash1);
	HashType bwHash2 = hash(revhash2);
	HashType bwHash3 = hash(revhash3);
	HashType bwHash4 = hash(revhash4);
	auto bwfound1 = list.hashToNode.find(bwHash1);
	auto bwfound2 = list.hashToNode.find(bwHash2);
	auto bwfound3 = list.hashToNode.find(bwHash3);
	auto bwfound4 = list.hashToNode.find(bwHash4);
	auto foundReal = found1;
	if (found1 != list.hashToNode.end())
	{
		list.addHashCharacterLength(sequenceCharacterLength, found1->second.second, seqCharLenStart, seqCharLenEnd, found1->second.first);
		foundReal = found1;
	}
	else if (found2 != list.hashToNode.end())
	{
		list.addHashCharacterLength(sequenceCharacterLength, found2->second.second, seqCharLenStart, seqCharLenEnd, found2->second.first);
		foundReal = found2;
	}
	else if (found3 != list.hashToNode.end())
	{
		list.addHashCharacterLength(sequenceCharacterLength, found3->second.second, seqCharLenStart, seqCharLenEnd, found3->second.first);
		foundReal = found3;
	}
	else if (found4 != list.hashToNode.end())
	{
		list.addHashCharacterLength(sequenceCharacterLength, found4->second.second, seqCharLenStart, seqCharLenEnd, found4->second.first);
		foundReal = found4;
	}
	if (foundReal != list.hashToNode.end())
	{
		auto found = foundReal->second;
		if (found1 == list.hashToNode.end()) list.hashToNode[fwHash1] = std::make_pair(found.first, found.second);
		if (found2 == list.hashToNode.end()) list.hashToNode[fwHash2] = std::make_pair(found.first, found.second);
		if (found3 == list.hashToNode.end()) list.hashToNode[fwHash3] = std::make_pair(found.first, found.second);
		if (found4 == list.hashToNode.end()) list.hashToNode[fwHash4] = std::make_pair(found.first, found.second);
		if (bwfound1 == list.hashToNode.end()) list.hashToNode[bwHash1] = std::make_pair(found.first, !found.second);
		if (bwfound2 == list.hashToNode.end()) list.hashToNode[bwHash2] = std::make_pair(found.first, !found.second);
		if (bwfound3 == list.hashToNode.end()) list.hashToNode[bwHash3] = std::make_pair(found.first, !found.second);
		if (bwfound4 == list.hashToNode.end()) list.hashToNode[bwHash4] = std::make_pair(found.first, !found.second);
		return std::make_pair(found, fwHash1);
	}

	assert(foundReal == list.hashToNode.end());
	assert(bwfound1 == list.hashToNode.end());
	assert(bwfound2 == list.hashToNode.end());
	assert(bwfound3 == list.hashToNode.end());
	assert(bwfound4 == list.hashToNode.end());
	size_t fwNode = list.size();
	list.hashToNode[fwHash1] = std::make_pair(fwNode, true);
	list.hashToNode[fwHash2] = std::make_pair(fwNode, true);
	list.hashToNode[fwHash3] = std::make_pair(fwNode, true);
	list.hashToNode[fwHash4] = std::make_pair(fwNode, true);
	list.hashToNode[bwHash1] = std::make_pair(fwNode, false);
	list.hashToNode[bwHash2] = std::make_pair(fwNode, false);
	list.hashToNode[bwHash3] = std::make_pair(fwNode, false);
	list.hashToNode[bwHash4] = std::make_pair(fwNode, false);
	list.addHashSequenceRLE(realSequence, fwHash1, previousHash, overlap);
	list.addHashCharacterLength(sequenceCharacterLength, seqCharLenStart, seqCharLenEnd, fwHash1, previousHash, overlap);
	assert(list.coverage.size() == fwNode);
	list.coverage.emplace_back(0);
	assert(list.edgeCoverage.size() == fwNode);
	list.edgeCoverage.emplace_back();
	assert(list.sequenceOverlap.size() == fwNode);
	list.sequenceOverlap.emplace_back();
	assert(list.fakeFwHashes.size() == fwNode);
	list.fakeFwHashes.emplace_back(fakeFwHash);
	assert(list.fakeBwHashes.size() == fwNode);
	list.fakeBwHashes.emplace_back(fakeBwHash);
	return std::make_pair(std::make_pair(fwNode, true), fwHash1);
}

template <typename F>
void findCollectedMinimizerPositions(const std::string& sequence, size_t kmerSize, const std::unordered_set<uint64_t>& collected, F callback)
{
	std::unordered_set<uint64_t> result;
	if (sequence.size() < kmerSize) return;
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	if (collected.count(fwkmerHasher.hash()) == 1)
	{
		callback(0, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
	}
	for (size_t i = kmerSize; i < sequence.size(); i++)
	{
		fwkmerHasher.addChar(sequence[i]);
		fwkmerHasher.removeChar(sequence[i-kmerSize]);
		if (collected.count(fwkmerHasher.hash()) == 1)
		{
			callback(i - kmerSize + 1, fwkmerHasher.getFwHash(), fwkmerHasher.getBwHash());
		}
	}
}

std::unordered_set<uint64_t> collectMinimizerHashes(const std::string& sequence, size_t kmerSize, size_t windowSize)
{
	std::unordered_set<uint64_t> result;
	if (sequence.size() < kmerSize + windowSize) return result;
	FastHasher fwkmerHasher { kmerSize };
	for (size_t i = 0; i < kmerSize; i++)
	{
		fwkmerHasher.addChar(sequence[i]);
	}
	std::deque<std::tuple<size_t, uint64_t>> minimizerOrder;
	minimizerOrder.emplace_back(0, fwkmerHasher.hash());
	for (size_t i = 0; i < windowSize-1; i++)
	{
		size_t seqPos = kmerSize+i;
		fwkmerHasher.addChar(sequence[seqPos]);
		fwkmerHasher.removeChar(sequence[seqPos-kmerSize]);
		size_t hash = fwkmerHasher.hash();
		while (minimizerOrder.size() > 0 && std::get<1>(minimizerOrder.back()) > hash) minimizerOrder.pop_back();
		minimizerOrder.emplace_back(i+1, hash);
	}
	auto pos = minimizerOrder.begin();
	while (pos != minimizerOrder.end() && std::get<1>(*pos) == std::get<1>(minimizerOrder.front()))
	{
		result.insert(std::get<1>(*pos));
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
				result.insert(std::get<1>(*pos));
				++pos;
			}
		}
		if (minimizerOrder.size() == 0 || hash == std::get<1>(minimizerOrder.front())) result.insert(fwkmerHasher.hash());
		minimizerOrder.emplace_back(i+1, hash);
	}
	return result;
}

HashList loadReadsAsPairedHashes(const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const bool hpc, const bool collapseRunLengths, const size_t pairingMinDistance, const size_t pairingMaxDistance)
{
	const size_t pairingMidDistance = (pairingMinDistance + pairingMaxDistance) / 2;
	std::unordered_set<uint64_t> collectedMinimizerHashes;
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&collectedMinimizerHashes, kmerSize, windowSize, hpc](const FastQ& read){
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
			auto collected = collectMinimizerHashes(seq, kmerSize, windowSize);
			collectedMinimizerHashes.insert(collected.begin(), collected.end());
		});
	}
	HashList paired { kmerSize, collapseRunLengths };
	size_t totalNodes = 0;
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&paired, &totalNodes, &collectedMinimizerHashes, kmerSize, windowSize, hpc, pairingMinDistance, pairingMaxDistance, pairingMidDistance](const FastQ& read){
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
			std::pair<size_t, bool> last { std::numeric_limits<size_t>::max(), true };
			HashType lastHash = 0;
			std::vector<size_t> minimizerPositions;
			std::vector<size_t> minimizerHash;
			findCollectedMinimizerPositions(seq, kmerSize, collectedMinimizerHashes, [&minimizerPositions, &minimizerHash](size_t pos, uint64_t fwHash, uint64_t bwHash)
			{
				minimizerPositions.push_back(pos);
				minimizerHash.push_back(std::min(fwHash, bwHash));
			});
			assert(minimizerPositions.size() >= 1);
			size_t previousMin = 0;
			size_t previousMid = 0;
			size_t previousMax = 0;
			size_t nextMin = 2;
			size_t nextMid = 2;
			size_t nextMax = 2;
			for (size_t i = 1; i < minimizerPositions.size()-1; i++)
			{
				if (minimizerPositions[i] - minimizerPositions[0] <= pairingMaxDistance + kmerSize) continue;
				if (minimizerPositions.back() - minimizerPositions[i] <= pairingMaxDistance + kmerSize) break;
				while (minimizerPositions[i] - minimizerPositions[previousMin + 1] > pairingMaxDistance + kmerSize) previousMin += 1;
				// if (minimizerPositions[i] - minimizerPositions[previousMin] <= pairingMaxDistance + kmerSize) continue;
				if (previousMid <= previousMin) previousMid = previousMin+1;
				while (minimizerPositions[i] - minimizerPositions[previousMid + 1] > pairingMidDistance + kmerSize) previousMid += 1;
				// if (minimizerPositions[i] - minimizerPositions[previousMid] <= pairingMidDistance + kmerSize) continue;
				if (previousMax <= previousMid) previousMax = previousMid+1;
				while (minimizerPositions[i] - minimizerPositions[previousMax + 1] > pairingMinDistance + kmerSize) previousMax += 1;
				// if (minimizerPositions[i] - minimizerPositions[previousMax] <= pairingMaxDistance + kmerSize) continue;
				assert(previousMin < i);
				assert(previousMid < i);
				assert(previousMax < i);
				assert(previousMin < previousMid);
				assert(previousMid < previousMax);
				assert(minimizerPositions[i] - minimizerPositions[previousMin] > pairingMaxDistance + kmerSize);
				assert(minimizerPositions[i] - minimizerPositions[previousMin+1] <= pairingMaxDistance + kmerSize);
				assert(minimizerPositions[i] - minimizerPositions[previousMid] > pairingMidDistance + kmerSize);
				assert(minimizerPositions[i] - minimizerPositions[previousMid+1] <= pairingMidDistance + kmerSize);
				assert(minimizerPositions[i] - minimizerPositions[previousMax] > pairingMinDistance + kmerSize);
				assert(minimizerPositions[i] - minimizerPositions[previousMax+1] <= pairingMinDistance + kmerSize);

				size_t previous = previousMin + 1;
				size_t previous2 = previousMid + 1;
				for (size_t j = previousMin+1; j <= previousMid; j++)
				{
					if (minimizerHash[j] < minimizerHash[previous]) previous = j;
				}
				for (size_t j = previousMid+1; j <= previousMax; j++)
				{
					if (minimizerHash[j] < minimizerHash[previous2]) previous2 = j;
				}
				assert(previous != previous2);

				if (nextMin <= i) nextMin = i+1;
				while (nextMin < minimizerPositions.size() && minimizerPositions[nextMin] - minimizerPositions[i] < pairingMinDistance + kmerSize) nextMin += 1;
				assert(nextMin < minimizerPositions.size());
				// if (nextMin == minimizerPositions.size()) break;
				if (nextMid <= nextMin) nextMid = nextMin+1;
				while (nextMid < minimizerPositions.size() && minimizerPositions[nextMid] - minimizerPositions[i] < pairingMidDistance + kmerSize) nextMid += 1;
				assert(nextMid < minimizerPositions.size());
				// if (nextMid == minimizerPositions.size()) break;
				if (nextMax <= nextMid) nextMax = nextMid+1;
				while (nextMax < minimizerPositions.size() && minimizerPositions[nextMax] - minimizerPositions[i] < pairingMaxDistance + kmerSize) nextMax += 1;
				// if (nextMax == minimizerPositions.size()) break;
				assert(nextMax < minimizerPositions.size());
				assert(nextMin < nextMid);
				assert(nextMid < nextMax);
				assert(nextMin > i);
				assert(minimizerPositions[nextMin] - minimizerPositions[i] >= pairingMinDistance + kmerSize);
				assert(minimizerPositions[nextMin-1] - minimizerPositions[i] < pairingMinDistance + kmerSize);
				assert(minimizerPositions[nextMid] - minimizerPositions[i] >= pairingMidDistance + kmerSize);
				assert(minimizerPositions[nextMid-1] - minimizerPositions[i] < pairingMidDistance + kmerSize);
				assert(minimizerPositions[nextMax] - minimizerPositions[i] >= pairingMaxDistance + kmerSize);
				assert(minimizerPositions[nextMax-1] - minimizerPositions[i] < pairingMaxDistance + kmerSize);

				size_t next = nextMin;
				size_t next2 = nextMid;
				for (size_t j = nextMin; j < nextMid; j++)
				{
					if (minimizerHash[j] < minimizerHash[next]) next = j;
				}
				for (size_t j = nextMid; j < nextMax; j++)
				{
					if (minimizerHash[j] < minimizerHash[next2]) next2 = j;
				}
				assert(next != next2);
				std::string hashseq1 = seq.substr(minimizerPositions[previous], kmerSize) + seq.substr(minimizerPositions[i], kmerSize) + seq.substr(minimizerPositions[next], kmerSize);
				std::string hashseq2 = seq.substr(minimizerPositions[previous2], kmerSize) + seq.substr(minimizerPositions[i], kmerSize) + seq.substr(minimizerPositions[next], kmerSize);
				std::string hashseq3 = seq.substr(minimizerPositions[previous], kmerSize) + seq.substr(minimizerPositions[i], kmerSize) + seq.substr(minimizerPositions[next2], kmerSize);
				std::string hashseq4 = seq.substr(minimizerPositions[previous2], kmerSize) + seq.substr(minimizerPositions[i], kmerSize) + seq.substr(minimizerPositions[next2], kmerSize);
				std::string revSeq1 = revCompRLE(hashseq1);
				std::string revSeq2 = revCompRLE(hashseq2);
				std::string revSeq3 = revCompRLE(hashseq3);
				std::string revSeq4 = revCompRLE(hashseq4);
				assert(hashseq1.size() == kmerSize * 3);
				assert(minimizerPositions[previous] + kmerSize < minimizerPositions[i]);
				assert(minimizerPositions[i] + kmerSize < minimizerPositions[next]);
				std::string_view hashview1 { hashseq1.data(), hashseq1.size() };
				std::string_view hashview2 { hashseq2.data(), hashseq2.size() };
				std::string_view hashview3 { hashseq3.data(), hashseq3.size() };
				std::string_view hashview4 { hashseq4.data(), hashseq4.size() };
				std::string_view realview { hashseq1.data() + kmerSize, kmerSize };
				std::string_view revhashview1 { revSeq1.data(), revSeq1.size() };
				std::string_view revhashview2 { revSeq2.data(), revSeq2.size() };
				std::string_view revhashview3 { revSeq3.data(), revSeq3.size() };
				std::string_view revhashview4 { revSeq4.data(), revSeq4.size() };
				std::string_view revrealview { revSeq1.data() + kmerSize, kmerSize };

				assert(last.first == std::numeric_limits<size_t>::max() || minimizerPositions[i] - minimizerPositions[i-1] <= kmerSize);
				std::pair<size_t, bool> current;
				size_t overlap = minimizerPositions[i-1] + kmerSize - minimizerPositions[i];
				std::tie(current, lastHash) = addNode(paired, hashseq1, hashseq2, hashseq3, hashseq4, realview, revhashview1, revhashview2, revhashview3, revhashview4, revrealview, lens, minimizerPositions[i], minimizerPositions[i] + kmerSize, lastHash, overlap, 0, 0);
				if (last.first != std::numeric_limits<size_t>::max() && minimizerPositions[i] - minimizerPositions[i-1] < kmerSize)
				{
					assert(minimizerPositions[i-1] + kmerSize >= minimizerPositions[i]);
					paired.addSequenceOverlap(last, current, overlap);
					auto pair = canon(last, current);
					paired.edgeCoverage[pair.first][pair.second] += 1;
				}
				paired.coverage[current.first] += 1;
				last = current;
				totalNodes += 1;
			}
		});
	}
	paired.buildReverseCompHashSequences();
	std::cerr << totalNodes << " nodes" << std::endl;
	std::cerr << paired.size() << " distinct fw/bw sequence nodes" << std::endl;
	return paired;
}

HashList loadReadsAsHashes(const std::vector<std::string>& files, const size_t kmerSize, const size_t windowSize, const bool hpc, const bool collapseRunLengths)
{
	HashList result { kmerSize, collapseRunLengths };
	size_t totalNodes = 0;
	std::unordered_set<uint64_t> collectedMinimizerHashes;
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&collectedMinimizerHashes, kmerSize, windowSize, hpc](const FastQ& read){
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
			auto collected = collectMinimizerHashes(seq, kmerSize, windowSize);
			collectedMinimizerHashes.insert(collected.begin(), collected.end());
		});
	}
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&result, &totalNodes, &collectedMinimizerHashes, kmerSize, windowSize, hpc](const FastQ& read){
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
			findCollectedMinimizerPositions(seq, kmerSize, collectedMinimizerHashes, [kmerSize, windowSize, &lastHash, &last, &lens, &seq, &revSeq, &result, &lastMinimizerPosition, &totalNodes](size_t pos, uint64_t fwHash, uint64_t bwHash)
			{
				assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos > lastMinimizerPosition);
				assert(lastMinimizerPosition == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= windowSize);
				assert(last.first == std::numeric_limits<size_t>::max() || pos - lastMinimizerPosition <= kmerSize);
				std::string_view minimizerSequence { seq.data() + pos, kmerSize };
				size_t revPos = seq.size() - (pos + kmerSize);
				std::string_view revMinimizerSequence { revSeq.data() + revPos, kmerSize };
				std::pair<size_t, bool> current;
				size_t overlap = lastMinimizerPosition + kmerSize - pos;
				std::tie(current, lastHash) = addNode(result, minimizerSequence, minimizerSequence, revMinimizerSequence, revMinimizerSequence, lens, pos, pos + kmerSize, lastHash, overlap, fwHash, bwHash);
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

void runMBG(const std::vector<std::string>& inputReads, const std::string& outputGraph, const size_t kmerSize, const size_t windowSize, const size_t minCoverage, const double minUnitigCoverage, const bool hpc, const bool collapseRunLengths, const size_t pairingMinDistance, const size_t pairingMaxDistance)
{
	auto beforeReading = getTime();
	HashList reads = (pairingMaxDistance == 0 ? loadReadsAsHashes(inputReads, kmerSize, windowSize, hpc, collapseRunLengths) : loadReadsAsPairedHashes(inputReads, kmerSize, windowSize, hpc, collapseRunLengths, pairingMinDistance, pairingMaxDistance));
	auto beforeUnitigs = getTime();
	auto unitigs = getUnitigGraph(reads, minCoverage);
	auto beforeFilter = getTime();
	if (minUnitigCoverage > minCoverage) unitigs = getUnitigs(filterUnitigsByCoverage(unitigs, minUnitigCoverage));
	auto beforeWrite = getTime();
	writeGraph(unitigs, outputGraph, reads);
	auto beforeStats = getTime();
	auto unitigStats = getSizeAndN50(reads, unitigs);
	auto afterStats = getTime();
	std::cerr << "reading and hashing sequences took " << formatTime(beforeReading, beforeUnitigs) << std::endl;
	std::cerr << "unitigifying took " << formatTime(beforeUnitigs, beforeFilter) << std::endl;
	std::cerr << "filtering unitigs took " << formatTime(beforeFilter, beforeWrite) << std::endl;
	std::cerr << "writing the graph took " << formatTime(beforeWrite, beforeStats) << std::endl;
	std::cerr << "calculating stats took " << formatTime(beforeStats, afterStats) << std::endl;
	std::cerr << "nodes: " << unitigs.numNodes() << std::endl;
	std::cerr << "edges: " << unitigs.numEdges() << std::endl;
	std::cerr << "assembly size " << unitigStats.first << " bp, N50 " << unitigStats.second << std::endl;
	std::cerr << "approximate number of k-mers ~ " << (unitigStats.first - unitigs.numNodes() * (kmerSize - windowSize/2 - 1)) << std::endl;
}
