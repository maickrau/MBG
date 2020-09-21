#include <cassert>
#include "LazyString.h"

LazyString::LazyString(TwobitView first, TwobitView second, size_t overlap) :
first(first),
second(second),
overlap(overlap)
{
	assert(overlap < first.size());
	assert(overlap < second.size());
}

char LazyString::operator[](size_t index) const
{
	if (index < first.size()) return first[index];
	assert(index >= first.size());
	size_t secondIndex = index - first.size() + overlap;
	assert(secondIndex < second.size());
	return second[secondIndex];
}

std::string_view LazyString::view(size_t start, size_t kmerSize)
{
	if (str.size() == 0)
	{
		str.reserve(size());
		str += first.toString();
		str += second.toSubstring(overlap);
	}
	assert(str.size() == size());
	return std::string_view { str.data() + start, kmerSize };
}

size_t LazyString::size() const
{
	return first.size() + second.size() - overlap;
}
