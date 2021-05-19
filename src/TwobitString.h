#ifndef TwobitString_h
#define TwobitString_h

#include <string>
#include <string_view>
#include <vector>
#include "VectorView.h"

class TwobitString
{
public:
	TwobitString();
	TwobitString(VectorView<uint16_t> str);
	uint16_t get(size_t i) const;
	void set(size_t i, uint16_t c);
	size_t size() const;
	TwobitString& operator=(const std::vector<uint16_t>& str);
	TwobitString& operator=(VectorView<uint16_t> str);
	std::vector<uint16_t> toString() const;
	void resize(size_t size);
	void push_back(uint16_t c);
	template <typename Iter>
	void insert(Iter start, Iter end)
	{
		Iter i = start;
		while (i != end)
		{
			push_back(*i);
			++i;
		}
	}
private:
	std::vector<uint16_t> data;
};


class TwobitView
{
public:
	TwobitView(const TwobitString& str, const size_t start, const size_t end);
	uint16_t operator[](size_t i) const;
	size_t size() const;
	std::vector<uint16_t> toString() const;
	std::vector<uint16_t> toSubstring(size_t substrStart) const;
private:
	const TwobitString& str;
	const size_t start;
	const size_t end;
};

TwobitString revCompMultiRLE(const TwobitString& str);

#endif
