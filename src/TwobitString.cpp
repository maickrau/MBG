#include <cassert>
#include "TwobitString.h"

TwobitString::TwobitString() :
data(),
realSize(0)
{}

TwobitString::TwobitString(const std::string_view& str) :
data(),
realSize(0)
{
	*this = str;
}

char TwobitString::get(size_t i) const
{
	size_t pos = i / 4;
	size_t off = (i % 4) * 2;
	return ((data[pos] >> off) & 3) + 1;
}

void TwobitString::set(size_t i, char c)
{
	assert(c >= 1);
	assert(c <= 4);
	size_t pos = i / 4;
	size_t off = (i % 4) * 2;
	data[pos] &= ~(3 << off);
	data[pos] |= (c-1) << off;
}

size_t TwobitString::size() const
{
	return realSize;
}

TwobitString& TwobitString::operator=(const std::string& str)
{
	data.clear();
	resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		set(i, str[i]);
	}
	return *this;
}

TwobitString& TwobitString::operator=(const std::string_view& str)
{
	data.clear();
	resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		set(i, str[i]);
	}
	return *this;
}

void TwobitString::resize(size_t size)
{
	size_t max = (size + 3) / 4;
	data.resize(max, 0);
	realSize = size;
}

void TwobitString::push_back(char c)
{
	size_t pos = realSize / 4;
	if (pos >= data.size())
	{
		data.resize(data.size() * 2, 0);
	}
	set(realSize, c);
	realSize += 1;
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
