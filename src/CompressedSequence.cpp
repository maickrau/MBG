#include <cassert>
#include <algorithm>
#include "CompressedSequence.h"
#include "ErrorMaskHelper.h"

CompressedSequence::CompressedSequence(const std::vector<uint16_t>& compressed, const std::vector<std::string>& expanded) :
	compressed(compressed),
	expanded(expanded)
{
	assert(compressed.size() == expanded.size());
}

uint16_t CompressedSequence::getCompressed(size_t i) const
{
	return compressed[i];
}

std::string CompressedSequence::getExpanded(size_t i) const
{
	return expanded[i];
}

size_t CompressedSequence::compressedSize() const
{
	return compressed.size();
}

void CompressedSequence::setCompressed(size_t i, uint16_t c)
{
	compressed[i] = c;
}

void CompressedSequence::setExpanded(size_t i, std::string seq)
{
	expanded[i] = seq;
}

std::vector<size_t> CompressedSequence::getExpandedPositions() const
{
	std::vector<size_t> result;
	result.resize(expanded.size()+1);
	result[0] = 0;
	for (size_t i = 0; i < expanded.size(); i++)
	{
		result[i+1] = result[i] + expanded[i].size();
	}
	return result;
}

std::string CompressedSequence::getExpandedSequence() const
{
	std::string result;
	for (size_t i = 0; i < expanded.size(); i++)
	{
		result += expanded[i];
	}
	return result;
}

CompressedSequence CompressedSequence::substr(size_t start, size_t len) const
{
	assert(start + len <= compressedSize());
	CompressedSequence result;
	result.compressed.insert(result.compressed.end(), compressed.begin() + start, compressed.begin() + start + len);
	result.expanded.insert(result.expanded.end(), expanded.begin() + start, expanded.begin() + start + len);
	return result;
}

CompressedSequence CompressedSequence::revComp() const
{
	CompressedSequence result = *this;
	result.compressed = revCompRLE(result.compressed);
	std::reverse(result.expanded.begin(), result.expanded.end());
	for (size_t i = 0; i < result.expanded.size(); i++)
	{
		result.expanded[i] = revCompRaw(result.expanded[i]);
	}
	return result;
}

void CompressedSequence::insertEnd(const CompressedSequence& seq)
{
	compressed.insert(compressed.end(), seq.compressed.begin(), seq.compressed.end());
	expanded.insert(expanded.end(), seq.expanded.begin(), seq.expanded.end());
}

void CompressedSequence::resize(size_t size)
{
	compressed.resize(size);
	expanded.resize(size);
}

std::vector<uint16_t> CompressedSequence::compressedSubstr(size_t start, size_t len) const
{
	assert(start + len <= compressed.size());
	return std::vector<uint16_t> { compressed.begin() + start, compressed.begin() + start + len };
}
