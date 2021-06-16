#ifndef ErrorMaskHelper_h
#define ErrorMaskHelper_h

#include <vector>
#include <string>
#include <tuple>

std::vector<std::pair<std::vector<uint16_t>, std::vector<uint8_t>>> multiRLECompress(const std::vector<uint16_t>& str, const size_t maxMaskLength);
size_t maxCode();
size_t getExpandedLength(uint16_t code, uint8_t length);
std::string getExpandedSequence(const std::vector<uint16_t>& seq, const std::vector<uint8_t>& length);
std::vector<uint16_t> revCompRLE(const std::vector<uint16_t>& str);
uint16_t complement(const uint16_t original);
size_t getOverlapFromRLE(const std::vector<std::pair<std::vector<uint16_t>, std::vector<uint8_t>>>& unitigSequences, std::pair<size_t, bool> fromUnitig, size_t rleOverlap);
std::vector<size_t> getRLEExpandedPositions(const std::vector<uint16_t>& seq, const std::vector<uint8_t>& lens);
std::pair<uint16_t, uint8_t> getCodeAndRunlength(const std::vector<uint16_t>& str, size_t start, size_t end, uint16_t motifLength);

#endif
