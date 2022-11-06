#include <cassert>
#include <cstring>
#include "MsatValueVector.h"

MsatValueVector::MsatValueChunk::MsatValueChunk() :
	vec(0),
	realsize(0),
	capacity(2)
{
}

MsatValueVector::MsatValueChunk::MsatValueChunk(MsatValueChunk&& other) :
	vec(other.vec),
	realsize(other.realsize),
	capacity(other.capacity)
{
	other.vec = 0;
	other.realsize = 0;
	other.capacity = 2;
}

MsatValueVector::MsatValueChunk::MsatValueChunk(const MsatValueChunk& other) :
	vec(0),
	realsize(0),
	capacity(2)
{
	*this = other;
}

MsatValueVector::MsatValueChunk::~MsatValueChunk()
{
	if (capacity >= 3) delete [] vec;
}

MsatValueVector::MsatValueChunk& MsatValueVector::MsatValueChunk::operator=(MsatValueChunk&& other)
{
	if (capacity >= 3) delete [] vec;
	vec = other.vec;
	realsize = other.realsize;
	capacity = other.capacity;
	other.vec = 0;
	other.realsize = 0;
	other.capacity = 2;
	return *this;
}

MsatValueVector::MsatValueChunk& MsatValueVector::MsatValueChunk::operator=(const MsatValueChunk& other)
{
	if (capacity >= 3) delete [] vec;
	if (other.capacity <= 2)
	{
		vec = other.vec;
		realsize = other.realsize;
		capacity = other.capacity;
		return *this;
	}
	vec = new uint32_t[other.capacity];
	realsize = other.realsize;
	capacity = other.capacity;
	memcpy(vec, other.vec, 4*realsize);
	return *this;
}

void MsatValueVector::MsatValueChunk::set(uint8_t index, uint16_t val)
{
	uint16_t got = get(index);
	if (got == val) return;
	if (got != 65535) erase(index);
	if (capacity <= 2 && realsize == 0)
	{
		vec = (uint32_t*)(((size_t)index << 24) + (size_t)val);
		realsize = 1;
		return;
	}
	if (capacity <= 2 && realsize == 1)
	{
		vec = (uint32_t*)(((size_t)vec & 0x00000000FFFFFFFF) + ((size_t)index << (size_t)56) + ((size_t)val << (size_t)32));
		realsize = 2;
		return;
	}
	if (capacity <= 2 && realsize == 2)
	{
		capacity = 5; // surprisingly 5 is empicirally lowest memory out of 4,5,6,7,8,9,10. you'd think otherwise but it is so.
		uint32_t* newVec = new uint32_t[capacity];
		realsize = 3;
		*(newVec+0) = (size_t)vec;
		*(newVec+1) = (size_t)vec >> 32;
		*(newVec+2) = ((uint32_t)index << 24) + (uint32_t)val;
		vec = newVec;
		return;
	}
	assert(capacity >= 3);
	assert(realsize < 255);
	if (realsize < capacity)
	{
		*(vec+realsize) = ((uint32_t)index << 24) + (uint32_t)val;
		realsize += 1;
		return;
	}
	assert(realsize == capacity);
	capacity *= 2;
	uint32_t* newVec = new uint32_t[capacity];
	memcpy(newVec, vec, 4*realsize);
	delete [] vec;
	vec = newVec;
	*(vec+realsize) = ((uint32_t)index << 24) + (uint32_t)val;
	realsize += 1;
}

void MsatValueVector::MsatValueChunk::erase(uint8_t index)
{
	uint16_t got = get(index);
	assert(got != 65535);
	assert(realsize >= 1);
	if (capacity <= 2 && realsize >= 1)
	{
		if ((((size_t)vec >> (size_t)24) & 255) == index)
		{
			vec = (uint32_t*)((size_t)vec >> (size_t)32);
			realsize -= 1;
			return;
		}
	}
	if (capacity <= 2 && realsize >= 2)
	{
		if ((((size_t)vec >> (size_t)56) & 255) == index)
		{
			realsize -= 1;
			return;
		}
	}
	assert(capacity >= 3);
	for (size_t i = 0; i < realsize; i++)
	{
		if ((((*(vec+i)) >> 24) & 255) == index)
		{
			for (size_t j = i+1; j < realsize; j++)
			{
				*(vec+j-1) = *(vec+j);
			}
			realsize -= 1;
			return;
		}
	}
	assert(false);
}

uint16_t MsatValueVector::MsatValueChunk::get(uint8_t index) const
{
	if (capacity <= 2)
	{
		if (realsize >= 1 && (((size_t)vec >> (size_t)24) & 255) == index)
		{
			return (size_t)vec & 65535;
		}
		if (realsize >= 2 && (((size_t)vec >> (size_t)56) & 255) == index)
		{
			return ((size_t)vec >> (size_t)32) & 65535;
		}
		return 65535;
	}
	for (size_t i = 0; i < realsize; i++)
	{
		if (((*(vec+i) >> 24) & 255) == index) return (*(vec+i)) & 65535;
	}
	return 65535;
}

size_t MsatValueVector::MsatValueChunk::size() const
{
	return realsize;
}

uint16_t MsatValueVector::get(size_t index) const
{
	size_t vecIndex = index / 256;
	size_t vecOffset = index % 256;
	return chunks[vecIndex].get(vecOffset);
}

void MsatValueVector::set(size_t index, uint16_t val)
{
	size_t vecIndex = index / 256;
	size_t vecOffset = index % 256;
	chunks[vecIndex].set(vecOffset, val);
}

void MsatValueVector::resize(size_t size)
{
	chunks.resize((size+255)/256);
}

void MsatValueVector::erase(size_t index)
{
	size_t vecIndex = index / 256;
	size_t vecOffset = index % 256;
	chunks[vecIndex].erase(vecOffset);
}
