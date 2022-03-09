#ifndef Serializer_h
#define Serializer_h

#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

namespace Serializer
{

	void write(std::ostream& stream, size_t value);
	void writeMostlyTwobits(std::ostream& stream, const std::vector<uint16_t>& value);
	template <typename T>
	void write(std::ostream& stream, const std::vector<T>& value)
	{
		write(stream, value.size());
		stream.write((const char*)value.data(), value.size() * sizeof(T));
	}
	void writeTwobits(std::ostream& stream, const std::string& value);
	void write(std::ostream& stream, const std::string& value);
	void read(std::istream& stream, size_t& value);
	void readMostlyTwobits(std::istream& stream, std::vector<uint16_t>& value);
	template <typename T>
	void read(std::istream& stream, std::vector<T>& value)
	{
		size_t size;
		read(stream, size);
		value.resize(size);
		stream.read((char*)value.data(), value.size() * sizeof(T));
	}
	void readTwobits(std::istream& stream, std::string& value);
	void read(std::istream& stream, std::string& value);
}

#endif
