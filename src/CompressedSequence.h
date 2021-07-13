#ifndef CompressedSequence_h
#define CompressedSequence_h

#include <vector>
#include <string>
#include <cstdint>

class CompressedSequence
{
public:
	using CompressedType = uint16_t;
	using ExpandedType = std::string;
	using CompressedIterator = typename std::vector<uint16_t>::const_iterator;
	using ExpandedIterator = typename std::vector<std::string>::const_iterator;
	void setCompressed(size_t i, uint16_t c);
	uint16_t getCompressed(size_t i) const;
	std::string getExpanded(size_t i) const;
	void setExpanded(size_t i, std::string seq);
	size_t compressedSize() const;
	CompressedIterator beginCompressed() const;
	CompressedIterator endCompressed() const;
	ExpandedIterator beginExpanded() const;
	ExpandedIterator endExpanded() const;
	std::vector<size_t> getExpandedPositions() const;
	std::string getExpandedSequence() const;
	CompressedSequence substr(size_t start, size_t len) const;
	CompressedSequence revComp() const;
	void insertEnd(CompressedIterator beginCompress, CompressedIterator endCompress, ExpandedIterator beginExpand, ExpandedIterator endExpand);
	void insertEnd(const CompressedSequence& seq);
	void resize(size_t size);
private:
	std::vector<uint16_t> compressed;
	std::vector<std::string> expanded;
};

#endif
