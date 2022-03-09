#include "Serializer.h"

namespace Serializer
{

	void write(std::ostream& stream, size_t value)
	{
		stream.write((const char*)&value, sizeof(size_t));
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

}
