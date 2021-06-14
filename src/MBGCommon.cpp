#include <string_view>
#include <algorithm>
#include <vector>
#include "MBGCommon.h"

HashType hash(std::string_view sequence)
{
	size_t half = sequence.size() / 2;
	size_t low = std::hash<std::string_view>{}(std::string_view { sequence.begin(), half });
	size_t high = std::hash<std::string_view>{}(std::string_view { (sequence.begin() + half), sequence.size() - half });
	return (HashType)low + (((HashType)high) << 64);
}

std::ostream& operator<<(std::ostream& os, HashType t)
{
	if (t == 0)
	{
		os << "0";
		return os;
	}
	std::string decimal;
	while (t != 0)
	{
		decimal += "0123456789"[t % 10];
		t /= 10;
	}
	std::reverse(decimal.begin(), decimal.end());
	os << decimal;
	return os;
}

std::istream& operator>>(std::istream& is, HashType& t)
{
	std::string decimal;
	is >> decimal;
	t = 0;
	for (size_t i = 0; i < decimal.size(); i++)
	{
		t *= 10;
		if (decimal[i] >= '0' && decimal[i] <= '9') t += decimal[i]-'0';
	}
	return is;
}

std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos)
{
	return std::make_pair(pos.first, !pos.second);
}

std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (to.first < from.first)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	if (to.first == from.first && !to.second && !from.second)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	return std::make_pair(from, to);
}

std::string revCompRLE(const std::string& original)
{
	std::string result { original.rbegin(), original.rend() };
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = complement(result[i]);
	}
	return result;
}

unsigned char complement(const unsigned char original)
{
	static std::vector<char> mapping { 0, 4, 3, 2, 1, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 25, 22, 19, 28, 21, 18, 27, 24, 17, 26, 23, 20 };
	return mapping[original];
}
