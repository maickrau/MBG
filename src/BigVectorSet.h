#ifndef BigVectorSet_h
#define BigVectorSet_h
#include <vector>
#include <cstddef>

class BigVectorSet
{
public:
	using const_iterator = std::vector<size_t>::const_iterator;
	bool get(size_t index) const;
	void set(size_t index);
	const_iterator begin() const;
	const_iterator end() const;
	void clear();
	void resize(size_t newSize);
	size_t size() const;
	size_t activeSize() const;
private:
	std::vector<bool> setbits;
	std::vector<size_t> active;
};

#endif
