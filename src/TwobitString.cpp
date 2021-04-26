#include <cassert>
#include "TwobitString.h"

TwobitString::TwobitString() :
data8bit()
// data(),
// realSize(0)
{}

TwobitString::TwobitString(const std::string_view& str) :
data8bit()
// data(),
// realSize(0)
{
	*this = str;
}

char TwobitString::get(size_t i) const
{
	return data8bit[i];
	// size_t pos = i / 4;
	// size_t off = (i % 4) * 2;
	// return ((data[pos] >> off) & 3) + 1;
}

void TwobitString::set(size_t i, char c)
{
	data8bit[i] = c;
	// assert(c >= 1);
	// assert(c <= 4);
	// size_t pos = i / 4;
	// size_t off = (i % 4) * 2;
	// data[pos] &= ~(3 << off);
	// data[pos] |= (c-1) << off;
}

size_t TwobitString::size() const
{
	return data8bit.size();
	// return realSize;
}

TwobitString& TwobitString::operator=(const std::string& str)
{
	data8bit.clear();
	data8bit = str;
	// data.clear();
	// resize(str.size());
	// for (size_t i = 0; i < str.size(); i++)
	// {
	// 	set(i, str[i]);
	// }
	return *this;
}

TwobitString& TwobitString::operator=(const std::string_view& str)
{
	data8bit.clear();
	data8bit = str;
	// data.clear();
	// resize(str.size());
	// for (size_t i = 0; i < str.size(); i++)
	// {
	// 	set(i, str[i]);
	// }
	return *this;
}

void TwobitString::resize(size_t size)
{
	data8bit.resize(size, 0);
	// size_t max = (size + 3) / 4;
	// data.resize(max, 0);
	// realSize = size;
}

void TwobitString::push_back(char c)
{
	data8bit.push_back(c);
	// size_t pos = realSize / 4;
	// if (pos >= data.size())
	// {
	// 	data.resize(data.size() * 2, 0);
	// }
	// set(realSize, c);
	// realSize += 1;
}

std::string TwobitString::toString() const
{
	std::string result;
	result.reserve(size());
	for (size_t i = 0; i < size(); i++)
	{
		result.push_back(get(i));
	}
	return result;
}

TwobitView::TwobitView(const TwobitString& str, const size_t start, const size_t end) :
str(str),
start(start),
end(end)
{}

char TwobitView::operator[](size_t i) const
{
	return str.get(start+i);
}

size_t TwobitView::size() const
{
	return end-start;
}

std::string TwobitView::toString() const
{
	std::string result;
	result.reserve(end - start);
	for (size_t i = start; i < end; i++)
	{
		result.push_back(str.get(i));
	}
	return result;
}

std::string TwobitView::toSubstring(size_t substrStart) const
{
	assert(start + substrStart < end);
	std::string result;
	result.reserve(end - (start + substrStart));
	for (size_t i = start + substrStart; i < end; i++)
	{
		result.push_back(str.get(i));
	}
	return result;
}

TwobitString revCompRLE(const TwobitView& original)
{
	static char mapping[5] { 0, 4, 3, 2, 1 };
	TwobitString result;
	result.resize(original.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result.set(i, mapping[(int)original[original.size()-1-i]]);
	}
	return result;
}

TwobitString revCompRLE(const TwobitString& original)
{
	static char mappingDinuc[29] { 0, 4, 3, 2, 1, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 25, 22, 19, 28, 21, 18, 27, 24, 17, 26, 23, 20 };
	// static char mapping[5] { 0, 4, 3, 2, 1 };
	TwobitString result;
	result.resize(original.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		assert((int)original.get(original.size()-1-i) >= 1);
		assert((int)original.get(original.size()-1-i) <= 28);
		result.set(i, mappingDinuc[(int)original.get(original.size()-1-i)]);
		// result.set(i, mapping[(int)original.get(original.size()-1-i)]);
	}
	return result;
}