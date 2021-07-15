#ifndef StringIndex_h
#define StringIndex_h

#include <vector>
#include <string>
#include <phmap.h>

class StringIndex
{
public:
	void init(size_t maxCode);
	uint32_t getIndex(uint16_t compressed, const std::string& expanded);
	std::string getString(uint16_t compressed, uint32_t index) const;
	void buildReverseIndex();
private:
	std::vector<phmap::flat_hash_map<std::string, uint32_t>> index;
	std::vector<std::vector<std::string>> reverseIndex;
};

#endif
