#ifndef MBGCommon_h
#define MBGCommon_h

#include <fstream>
#include <tuple>
#include <vector>
#include "VectorView.h"

using HashType = unsigned __int128;
using NodeType = size_t;

HashType hash(std::string_view sequence);
HashType hash(VectorView<uint16_t> sequence);
HashType hash(std::vector<uint16_t> sequence);
std::ostream& operator<<(std::ostream& os, HashType t);
std::istream& operator>>(std::istream& is, HashType& t);
std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos);
std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to);

namespace std
{
	template <> struct hash<HashType>
	{
		size_t operator()(HashType x) const
		{
			return (size_t)x ^ (size_t)(x >> 64);
		}
	};
	template <> struct hash<std::pair<HashType, bool>>
	{
		size_t operator()(std::pair<HashType, bool> x) const
		{
			return hash<HashType>{}(x.first);
		}
	};
	template <> struct hash<std::pair<HashType, HashType>>
	{
		size_t operator()(std::pair<HashType, HashType> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
	template <> struct hash<std::pair<size_t, bool>>
	{
		size_t operator()(std::pair<size_t, bool> x) const
		{
			return (size_t)x.first;
		}
	};
	template <> struct hash<std::pair<size_t, size_t>>
	{
		size_t operator()(std::pair<size_t, size_t> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
}

#endif
