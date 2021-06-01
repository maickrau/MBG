#ifndef MultiRLE_h
#define MultiRLE_h

#include <vector>
#include <string>
#include <tuple>

std::vector<std::pair<std::vector<uint16_t>, std::vector<uint16_t>>> multiRLECompress(const std::string& str);
std::string multiRLEDecompress(const std::vector<uint16_t>& codes, const std::vector<uint16_t>& lens);
size_t maxCode();
std::vector<uint16_t> revCompMultiRLE(const std::vector<uint16_t>& codes);
uint16_t reverseComplement(uint16_t code);
size_t getExpandedLength(uint16_t code, uint16_t length);

#endif
