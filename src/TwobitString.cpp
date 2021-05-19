#include <cassert>
#include "TwobitString.h"
#include "MultiRLE.h"

TwobitString::TwobitString() :
data()
{}

TwobitString::TwobitString(VectorView<uint16_t> str) :
data()
{
	*this = str;
}

uint16_t TwobitString::get(size_t i) const
{
	return data[i];
}

void TwobitString::set(size_t i, uint16_t c)
{
	data[i] = c;
}

size_t TwobitString::size() const
{
	return data.size();
}

TwobitString& TwobitString::operator=(const std::vector<uint16_t>& str)
{
	data.clear();
	resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		set(i, str[i]);
	}
	return *this;
}

TwobitString& TwobitString::operator=(VectorView<uint16_t> str)
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
	data.resize(size);
}

void TwobitString::push_back(uint16_t c)
{
	data.push_back(c);
}

std::vector<uint16_t> TwobitString::toString() const
{
	std::vector<uint16_t> result;
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

uint16_t TwobitView::operator[](size_t i) const
{
	return str.get(start+i);
}

size_t TwobitView::size() const
{
	return end-start;
}

std::vector<uint16_t> TwobitView::toString() const
{
	std::vector<uint16_t> result;
	result.reserve(end - start);
	for (size_t i = start; i < end; i++)
	{
		result.push_back(str.get(i));
	}
	return result;
}

std::vector<uint16_t> TwobitView::toSubstring(size_t substrStart) const
{
	assert(start + substrStart < end);
	std::vector<uint16_t> result;
	result.reserve(end - (start + substrStart));
	for (size_t i = start + substrStart; i < end; i++)
	{
		result.push_back(str.get(i));
	}
	return result;
}

TwobitString revCompMultiRLE(const TwobitString& str)
{
	TwobitString result;
	result.resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		result.set(i, reverseComplement(str.get(str.size()-i-1)));
	}
	return result;
}

