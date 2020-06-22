#include <cassert>
#include "CommonUtils.h"

namespace CommonUtils
{
	std::string ReverseComplement(std::string str)
	{
		std::string result;
		result.reserve(str.size());
		for (int i = str.size()-1; i >= 0; i--)
		{
			result += Complement(str[i]);
		}
		return result;
	}

	char Complement(char c)
	{
		switch (c)
		{
			case 'A':
			case 'a':
				return 'T';
			case 'C':
			case 'c':
				return 'G';
			case 'T':
			case 't':
				return 'A';
			case 'G':
			case 'g':
				return 'C';
			case 'N':
			case 'n':
				return 'N';
			case 'U':
			case 'u':
				return 'A';
			case 'R':
			case 'r':
				return 'Y';
			case 'Y':
			case 'y':
				return 'R';
			case 'K':
			case 'k':
				return 'M';
			case 'M':
			case 'm':
				return 'K';
			case 'S':
			case 's':
				return 'S';
			case 'W':
			case 'w':
				return 'W';
			case 'B':
			case 'b':
				return 'V';
			case 'V':
			case 'v':
				return 'B';
			case 'D':
			case 'd':
				return 'H';
			case 'H':
			case 'h':
				return 'D';
			default:
				assert(false);
				return 'N';
		}
	}

}
