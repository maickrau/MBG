#include <algorithm>
#include "CompressedSequence.h"
#include "ErrorMaskHelper.h"

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

CompressedSequence::CompressedIterator CompressedSequence::beginCompressed() const
{
	return compressed.begin();
}

CompressedSequence::CompressedIterator CompressedSequence::endCompressed() const
{
	return compressed.end();
}

CompressedSequence::ExpandedIterator CompressedSequence::beginExpanded() const
{
	return expanded.begin();
}

CompressedSequence::ExpandedIterator CompressedSequence::endExpanded() const
{
	return expanded.end();
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
	CompressedSequence result;
	result.insertEnd(beginCompressed() + start, beginCompressed() + start + len, beginExpanded() + start, beginExpanded() + start + len);
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

void CompressedSequence::insertEnd(CompressedIterator beginCompress, CompressedIterator endCompress, ExpandedIterator beginExpand, ExpandedIterator endExpand)
{
	compressed.insert(compressed.end(), beginCompress, endCompress);
	expanded.insert(expanded.end(), beginExpand, endExpand);
}

void CompressedSequence::insertEnd(const CompressedSequence& seq)
{
	insertEnd(seq.beginCompressed(), seq.endCompressed(), seq.beginExpanded(), seq.endExpanded());
}

void CompressedSequence::resize(size_t size)
{
	compressed.resize(size);
	expanded.resize(size);
}
