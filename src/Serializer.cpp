#include <iostream>
#include <cassert>
#include "Serializer.h"

void writeTightly(std::vector<uint8_t>& compressed, size_t pos, size_t numBits, size_t value)
{
	while (true)
	{
		size_t bitsInByte = 8 - (pos % 8);
		if (numBits <= bitsInByte)
		{
			assert((value << (pos % 8)) < 256);
			compressed[pos/8] |= (uint8_t)(value << (pos % 8));
			return;
		}
		size_t addHere = (value >> (numBits - bitsInByte));
		assert((addHere << (pos % 8)) < 256);
		compressed[pos/8] |= (uint8_t)(addHere << (pos % 8));
		value &= (1 << (numBits - bitsInByte)) - 1;
		numBits -= bitsInByte;
		pos += bitsInByte;
	}
}

void writeTightly(std::vector<uint8_t>& compressed, size_t pos, uint16_t value)
{
	if (value < 4)
	{
		writeTightly(compressed, pos+1, 2, value);
	}
	else
	{
		compressed[pos/8] |= (uint8_t)(1 << (pos % 8));
		writeTightly(compressed, pos+1, 16, value);
	}
}

uint16_t readTightly(std::vector<uint8_t>& compressed, size_t pos)
{
	bool longValue = ((compressed[pos/8] >> (pos % 8)) & 1) == 1;
	pos += 1;
	uint16_t value = 0;
	size_t numBits = longValue ? 16 : 2;
	while (numBits > 0)
	{
		size_t bitsInByte = 8 - (pos % 8);
		value <<= std::min(numBits, bitsInByte);
		if (numBits <= bitsInByte)
		{
			value |= ((size_t)compressed[pos/8] >> (pos % 8)) & ((1 << numBits) - 1);
			return value;
		}
		value |= ((size_t)compressed[pos/8] >> (pos % 8));
		pos += bitsInByte;
		numBits -= bitsInByte;
	}
	assert(false);
	return value;
}

namespace Serializer
{

	void write(std::ostream& stream, size_t value)
	{
		stream.write((const char*)&value, sizeof(size_t));
	}

	void writeTwobits(std::ostream& stream, const std::string& value)
	{
		write(stream, value.size());
		std::vector<uint8_t> twobits;
		twobits.resize((value.size() + 3) / 4, 0);
		for (size_t i = 0; i < value.size(); i++)
		{
			uint8_t val = 0;
			switch(value[i])
			{
				case 'a':
				case 'A':
					val = 0;
					break;
				case 'c':
				case 'C':
					val = 1;
					break;
				case 'g':
				case 'G':
					val = 2;
					break;
				case 't':
				case 'T':
					val = 3;
					break;
				default:
					assert(false);
					break;
			}
			twobits[i / 4] |= val << ((i % 4) * 2);
		}
		stream.write((const char*)twobits.data(), twobits.size());
	}

	void writeMonotoneIncreasing(std::ostream& stream, const std::vector<uint64_t>& value)
	{
		size_t numBytes = 0;
		for (size_t i = 0; i < value.size(); i++)
		{
			size_t val = value[i];
			if (i > 0)
			{
				assert(value[i] >= value[i-1]);
				val = value[i] - value[i-1];
			}
			if (val < 128) // 2^7, one byte minus one bit for length
			{
				numBytes += 1;
			}
			else if (val < 16384) // 2^14, two bytes minus two bits for length
			{
				numBytes += 2;
			}
			else
			{
				assert(val < std::numeric_limits<uint64_t>::max() / 4); // minus two bits for length
				numBytes += 8;
			}
		}
		std::vector<uint8_t> compressed;
		compressed.resize(numBytes, 0);
		write(stream, value.size());
		write(stream, compressed.size());
		size_t pos = 0;
		for (size_t i = 0; i < value.size(); i++)
		{
			size_t val = value[i];
			if (i > 0) val = value[i] - value[i-1];
			if (val < 128)
			{
				compressed[pos] = val;
				assert((compressed[pos] & 128) == 0);
				pos += 1;
			}
			else if (val < 16384)
			{
				assert((val >> 7) < 128);
				compressed[pos] = 128 + (val >> 7);
				compressed[pos+1] = val & 127;
				assert((compressed[pos] & 128) == 128);
				assert((compressed[pos+1] & 128) == 0);
				pos += 2;
			}
			else
			{
				assert((val >> 55) < 128);
				compressed[pos] = 128 + (val >> 55);
				compressed[pos+1] = 128 + ((val >> 48) & 127);
				compressed[pos+2] = (val >> 40) & 255;
				compressed[pos+3] = (val >> 32) & 255;
				compressed[pos+4] = (val >> 24) & 255;
				compressed[pos+5] = (val >> 16) & 255;
				compressed[pos+6] = (val >> 8) & 255;
				compressed[pos+7] = (val) & 255;
				assert((compressed[pos] & 128) == 128);
				assert((compressed[pos+1] & 128) == 128);
				pos += 8;
			}
		}
		assert(pos == numBytes);
		stream.write((const char*)compressed.data(), compressed.size());
	}

	void writeMostlyTwobits(std::ostream& stream, const std::vector<uint16_t>& value)
	{
		write(stream, value.size());
		std::vector<uint8_t> compressed;
		size_t numBits = 0;
		for (size_t i = 0; i < value.size(); i++)
		{
			if (value[i] < 4)
			{
				numBits += 3;
			}
			else
			{
				numBits += 17;
			}
		}
		compressed.resize((numBits+7)/8, 0);
		write(stream, compressed.size());
		size_t pos = 0;
		for (size_t i = 0; i < value.size(); i++)
		{
			writeTightly(compressed, pos, value[i]);
			if (value[i] < 4)
			{
				pos += 3;
			}
			else
			{
				pos += 17;
			}
		}
		assert(pos == numBits);
		stream.write((const char*)compressed.data(), compressed.size());
	}

	void write(std::ostream& stream, const std::string& value)
	{
		write(stream, value.size());
		stream.write((const char*)value.data(), value.size());
	}

	void read(std::istream& stream, size_t& value)
	{
		stream.read((char*)&value, sizeof(size_t));
	}

	void read(std::istream& stream, std::string& value)
	{
		size_t size;
		read(stream, size);
		value.resize(size);
		stream.read((char*)value.data(), value.size());
	}

	void readTwobits(std::istream& stream, std::string& value)
	{
		size_t size;
		read(stream, size);
		std::vector<uint8_t> twobits;
		twobits.resize((size+3)/4, 0);
		stream.read((char*)twobits.data(), twobits.size());
		value.resize(size);
		for (size_t i = 0; i < value.size(); i++)
		{
			uint8_t val = (twobits[i/4] >> ((i % 4) * 2)) & 3;
			switch(val)
			{
				case 0:
					value[i] = 'A';
					break;
				case 1:
					value[i] = 'C';
					break;
				case 2:
					value[i] = 'G';
					break;
				case 3:
					value[i] = 'T';
					break;
				default:
					assert(false);
					break;
			}
		}
	}

	void readMonotoneIncreasing(std::istream& stream, std::vector<uint64_t>& value)
	{
		size_t realSize;
		size_t numBytes;
		read(stream, realSize);
		read(stream, numBytes);
		value.resize(realSize);
		std::vector<uint8_t> compressed;
		compressed.resize(numBytes, 0);
		stream.read((char*)compressed.data(), numBytes);
		size_t pos = 0;
		size_t runningTotal = 0;
		for (size_t i = 0; i < value.size(); i++)
		{
			size_t val = 0;
			assert(pos < numBytes);
			if ((compressed[pos] & 128) == 0)
			{
				val = compressed[pos];
				pos += 1;
			}
			else if ((compressed[pos+1] & 128) == 0)
			{
				val = compressed[pos] & 127;
				val <<= 7;
				val += compressed[pos+1];
				pos += 2;
			}
			else
			{
				val = compressed[pos] & 127;
				val <<= 7;
				val += compressed[pos+1] & 127;
				val <<= 7;
				val += compressed[pos+2];
				val <<= 8;
				val += compressed[pos+3];
				val <<= 8;
				val += compressed[pos+4];
				val <<= 8;
				val += compressed[pos+5];
				val <<= 8;
				val += compressed[pos+6];
				val <<= 8;
				val += compressed[pos+7];
				pos += 8;
			}
			runningTotal += val;
			value[i] = runningTotal;
		}
		assert(pos == numBytes);
	}

	void readMostlyTwobits(std::istream& stream, std::vector<uint16_t>& value)
	{
		size_t realSize;
		size_t numBytes;
		read(stream, realSize);
		read(stream, numBytes);
		std::vector<uint8_t> compressed;
		compressed.resize(numBytes, 0);
		stream.read((char*)compressed.data(), compressed.size());
		size_t pos = 0;
		value.resize(realSize);
		for (size_t i = 0; i < value.size(); i++)
		{
			uint16_t val = readTightly(compressed, pos);
			if (val < 4)
			{
				pos += 3;
			}
			else
			{
				pos += 17;
			}
			value[i] = val;
		}
		assert((pos + 7) / 8 == numBytes);
	}

}
