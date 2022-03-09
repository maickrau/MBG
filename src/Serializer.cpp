#include <cassert>
#include "Serializer.h"

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

}
