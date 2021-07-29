#include <cassert>
#include <algorithm>
#include "CompressedSequence.h"
#include "ErrorMaskHelper.h"

CompressedSequence::CompressedSequence(const std::vector<uint16_t>& compressed, const std::vector<uint32_t>& expanded) :
	compressed(compressed)
{
	assert(compressed.size() == expanded.size());
	simpleExpanded.resize(expanded.size());
	for (size_t i = 0; i < expanded.size(); i++)
	{
		if (expanded[i] < 256)
		{
			simpleExpanded[i] = expanded[i];
		}
		else
		{
			complexExpanded[i] = expanded[i];
		}
	}
	for (size_t i = 0; i < expanded.size(); i++)
	{
		assert(getExpanded(i) > 0 || getCompressed(i) >= 4);
	}
}

uint16_t CompressedSequence::getCompressed(size_t i) const
{
	assert(i < compressed.size());
	return compressed[i];
}

uint32_t CompressedSequence::getExpanded(size_t i) const
{
	assert(i < simpleExpanded.size());
	auto found = complexExpanded.find(i);
	if (found != complexExpanded.end()) return found->second;
	return simpleExpanded[i];
}

std::string CompressedSequence::getExpandedStr(size_t i, const StringIndex& index) const
{
	return index.getString(getCompressed(i), getExpanded(i));
}

size_t CompressedSequence::compressedSize() const
{
	return compressed.size();
}

void CompressedSequence::setCompressed(size_t i, uint16_t c)
{
	assert(i < compressed.size());
	compressed[i] = c;
}

void CompressedSequence::setExpanded(size_t i, uint32_t seq)
{
	assert(i < simpleExpanded.size());
	if (seq < 256)
	{
		auto found = complexExpanded.find(i);
		if (found != complexExpanded.end()) complexExpanded.erase(found);
		simpleExpanded[i] = seq;
		return;
	}
	complexExpanded[i] = seq;
}

std::vector<size_t> CompressedSequence::getExpandedPositions(const StringIndex& index) const
{
	std::vector<size_t> result;
	result.resize(simpleExpanded.size()+1);
	result[0] = 0;
	for (size_t i = 0; i < simpleExpanded.size(); i++)
	{
		result[i+1] = result[i] + index.getString(compressed[i], getExpanded(i)).size();
	}
	return result;
}

std::string CompressedSequence::getExpandedSequence(const StringIndex& index) const
{
	std::string result;
	for (size_t i = 0; i < simpleExpanded.size(); i++)
	{
		result += index.getString(compressed[i], getExpanded(i));
	}
	return result;
}

CompressedSequence CompressedSequence::substr(size_t start, size_t len) const
{
	assert(start + len <= compressedSize());
	CompressedSequence result;
	result.compressed.insert(result.compressed.end(), compressed.begin() + start, compressed.begin() + start + len);
	result.simpleExpanded.insert(result.simpleExpanded.end(), simpleExpanded.begin() + start, simpleExpanded.begin() + start + len);
	for (auto pair : complexExpanded)
	{
		if (pair.first >= start && pair.first < start + len) result.complexExpanded[pair.first - start] = pair.second;
	}
	return result;
}

CompressedSequence CompressedSequence::revComp(const StringIndex& stringIndex) const
{
	CompressedSequence result;
	result.compressed.resize(compressed.size());
	result.simpleExpanded.resize(simpleExpanded.size());
	for (size_t i = 0; i < compressed.size(); i++)
	{
		uint16_t comp = complement(compressed[i]);
		uint32_t expanded = stringIndex.getReverseIndex(comp, getExpanded(i));
		result.setCompressed(compressed.size()-1-i, comp);
		result.setExpanded(compressed.size()-1-i, expanded);
	}
	return result;
}

void CompressedSequence::insertEnd(const CompressedSequence& seq)
{
	size_t oldSize = compressedSize();
	compressed.insert(compressed.end(), seq.compressed.begin(), seq.compressed.end());
	simpleExpanded.insert(simpleExpanded.end(), seq.simpleExpanded.begin(), seq.simpleExpanded.end());
	for (auto pair : seq.complexExpanded)
	{
		complexExpanded[oldSize + pair.first] = pair.second;
	}
}

void CompressedSequence::resize(size_t size)
{
	compressed.resize(size);
	simpleExpanded.resize(size);
}

std::vector<uint16_t> CompressedSequence::compressedSubstr(size_t start, size_t len) const
{
	assert(start + len <= compressed.size());
	return std::vector<uint16_t> { compressed.begin() + start, compressed.begin() + start + len };
}
