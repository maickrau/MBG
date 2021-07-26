#ifndef RankBitvector_h
#define RankBitvector_h

#include <cstdint>
#include <limits>
#include <vector>

class RankBitvector
{
private:
	static constexpr size_t BitsPerChunk = 64;
	static constexpr size_t SmallRanksPerBig = std::numeric_limits<uint16_t>::max() / 64;
public:
	RankBitvector(size_t size);
	void set(size_t i, bool value);
	bool get(size_t i) const;
	void buildRanks();
	size_t getRank(size_t i) const;
	size_t size() const;
private:
	bool ranksBuilt;
	std::vector<uint64_t> bits;
	std::vector<uint16_t> smallRanks;
	std::vector<size_t> bigRanks;
	size_t realSize;
};

#endif
