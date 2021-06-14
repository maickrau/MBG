#include "ReadHelper.h"

std::pair<std::string, std::vector<uint8_t>> dinucRunLengthEncode(const std::string& rleString, const std::vector<uint8_t>& runLengths)
{
	// how to preserve the run lengths of the nucleotides inside the run length of a dinucleotide?
	// no solution yet, so remove the run lengths of the nucleotides outside dinucs as well so everywhere is equally biased
	std::string result;
	std::vector<uint8_t> lens;
	size_t i = 0;
	bool boundaryStart = false;
	for (; i <= rleString.size() - 4; i++)
	{
		if (rleString[i] != rleString[i+2] || rleString[i+1] != rleString[i+3])
		{
			result.push_back(rleString[i]);
			lens.push_back(1);
			boundaryStart = false;
			continue;
		}
		size_t rleEnd = i+4;
		if (boundaryStart) i += 1;
		while (rleEnd < rleString.size() && rleString[rleEnd] == rleString[rleEnd-2]) rleEnd += 1;
		size_t rleLen = rleEnd - i;
		char firstChar = rleString[i];
		char secondChar = rleString[i+1];
		assert(firstChar != secondChar);
		assert(firstChar >= 1 && firstChar <= 4);
		assert(secondChar >= 1 && secondChar <= 4);
		boundaryStart = false;
		i = rleEnd - 1;
		char boundaryChar = 0;
		if (i <= rleString.size()-4 && rleString[i] == rleString[i+2] && rleString[i+1] == rleString[i+3])
		{
			// two adjacent dimers share one bp
			// split the shared bp apart
			// but add it after the dimer
			boundaryChar = rleString[i];
			boundaryStart = true;
			i -= 1;
			rleLen -= 1;
		}
		// odd length, type ATATA
		if (rleLen % 2 == 1)
		{
			if (firstChar == 1 && secondChar == 2) result.push_back(5); // ACA
			if (firstChar == 1 && secondChar == 3) result.push_back(6); // AGA
			if (firstChar == 1 && secondChar == 4) result.push_back(7); // ATA
			if (firstChar == 2 && secondChar == 1) result.push_back(8); // CAC
			if (firstChar == 2 && secondChar == 3) result.push_back(9); // CGC
			if (firstChar == 2 && secondChar == 4) result.push_back(10); // CTC
			if (firstChar == 3 && secondChar == 1) result.push_back(11); // GAG
			if (firstChar == 3 && secondChar == 2) result.push_back(12); // GCG
			if (firstChar == 3 && secondChar == 4) result.push_back(13); // GTG
			if (firstChar == 4 && secondChar == 1) result.push_back(14); // TAT
			if (firstChar == 4 && secondChar == 2) result.push_back(15); // TCT
			if (firstChar == 4 && secondChar == 3) result.push_back(16); // TGT
			lens.push_back((rleLen - 1) / 2);
		}
		// even length, type ATAT
		else
		{
			if (firstChar == 1 && secondChar == 2) result.push_back(17); // AC
			if (firstChar == 1 && secondChar == 3) result.push_back(18); // AG
			if (firstChar == 1 && secondChar == 4) result.push_back(19); // AT
			if (firstChar == 2 && secondChar == 1) result.push_back(20); // CA
			if (firstChar == 2 && secondChar == 3) result.push_back(21); // CG
			if (firstChar == 2 && secondChar == 4) result.push_back(22); // CT
			if (firstChar == 3 && secondChar == 1) result.push_back(23); // GA
			if (firstChar == 3 && secondChar == 2) result.push_back(24); // GC
			if (firstChar == 3 && secondChar == 4) result.push_back(25); // GT
			if (firstChar == 4 && secondChar == 1) result.push_back(26); // TA
			if (firstChar == 4 && secondChar == 2) result.push_back(27); // TC
			if (firstChar == 4 && secondChar == 3) result.push_back(28); // TG
			lens.push_back(rleLen / 2);
		}
		if (boundaryChar != 0)
		{
			result.push_back(boundaryChar);
			lens.push_back(1);
		}
		assert(result.size() == lens.size());
	}
	assert(result.size() == lens.size());
	for (; i < rleString.size(); i++)
	{
		result.push_back(rleString[i]);
		lens.push_back(1);
	}
	assert(result.size() == lens.size());
	return std::make_pair(result, lens);
}
