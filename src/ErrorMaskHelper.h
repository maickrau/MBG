#ifndef ErrorMaskHelper_h
#define ErrorMaskHelper_h

#include <vector>
#include <string>
#include <tuple>
#include "MBGCommon.h"

std::vector<std::pair<SequenceCharType, SequenceLengthType>> multiRLECompress(const SequenceCharType& str, const size_t maxMaskLength);
size_t maxCode();
size_t getExpandedLength(CharType code, LengthType length);
std::string getExpandedSequence(const SequenceCharType& seq, const SequenceLengthType& length);
SequenceCharType revCompRLE(const SequenceCharType& str);
CharType complement(const CharType original);
size_t getOverlapFromRLE(const std::vector<std::pair<SequenceCharType, SequenceLengthType>>& unitigSequences, std::pair<size_t, bool> fromUnitig, size_t rleOverlap);
std::vector<size_t> getRLEExpandedPositions(const SequenceCharType& seq, const SequenceLengthType& lens);
std::pair<CharType, LengthType> getCodeAndRunlength(const SequenceCharType& str, size_t start, size_t end, uint16_t motifLength);

#endif
