#include <cassert>
#include "BigVectorSet.h"

bool BigVectorSet::get(size_t index) const
{
	return setbits[index];
}

void BigVectorSet::set(size_t index)
{
	if (setbits[index]) return;
	setbits[index] = true;
	active.push_back(index);
}

BigVectorSet::const_iterator BigVectorSet::begin() const
{
	return active.begin();
}

BigVectorSet::const_iterator BigVectorSet::end() const
{
	return active.end();
}

void BigVectorSet::clear()
{
	for (auto index : active)
	{
		assert(setbits[index]);
		setbits[index] = false;
	}
	active.clear();
}

void BigVectorSet::resize(size_t newSize)
{
	setbits.resize(newSize, false);
}

size_t BigVectorSet::size() const
{
	return setbits.size();
}

size_t BigVectorSet::activeSize() const
{
	return active.size();
}
