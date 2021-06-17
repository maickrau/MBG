#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>
#include "ErrorMaskHelper.h"

constexpr size_t MaxMotifLength = 6;

uint16_t getNumBefore(uint16_t motifLength)
{
	size_t result = 0;
	for (uint16_t i = 1; i < motifLength; i++)
	{
		result += pow(4, i) * i;
	}
	assert(result < std::numeric_limits<uint16_t>::max());
	return result;
}

CharType getReverseComplement(CharType code)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	code -= getNumBefore(motifLength);
	assert(code < pow(4, motifLength) * motifLength);
	uint16_t motif = code / motifLength;
	assert(motif < pow(4, motifLength));
	uint16_t overhang = code % motifLength;
	assert(overhang < motifLength);
	if (overhang > 0)
	{
		// rotate to start with the overhang
		motif = (motif >> ((motifLength - overhang) * 2)) + ((motif & ((1 << ((motifLength - overhang) * 2)) - 1)) << ((overhang) * 2));
	}
	assert(motif < pow(4, motifLength));
	uint16_t newMotif = 0;
	for (uint16_t i = 0; i < motifLength; i++)
	{
		newMotif <<= 2;
		newMotif |= (~motif) & 3;
		motif >>= 2;
	}
	assert(newMotif < pow(4, motifLength));
	CharType newCode = getNumBefore(motifLength) + newMotif * motifLength + overhang;
	return newCode;
}

std::vector<CharType> getReverseComplements()
{
	SequenceCharType result;
	for (uint16_t motifLength = 1; motifLength <= MaxMotifLength; motifLength++)
	{
		for (uint16_t i = 0; i < getNumBefore(motifLength+1) - getNumBefore(motifLength); i++)
		{
			result.push_back(getReverseComplement(getNumBefore(motifLength) + i));
		}
	}
	assert(result.size() == getNumBefore(MaxMotifLength+1));
	assert(result.size() < std::numeric_limits<CharType>::max());
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[result[i]] == i);
	}
	return result;
}

std::vector<CharType> multiRLEReverseComplements = getReverseComplements();

CharType complement(CharType code)
{
	assert(code < multiRLEReverseComplements.size());
	return multiRLEReverseComplements[code];
}

size_t maxCode()
{
	return multiRLEReverseComplements.size();
}

std::pair<CharType, LengthType> getCodeAndRunlength(const SequenceCharType& str, size_t start, size_t end, uint16_t motifLength)
{
	assert(end > start);
	assert(end >= start + motifLength);
	uint16_t overhang = (end - start) % motifLength;
	uint16_t motif = 0;
	for (size_t i = start; i < start+motifLength; i++)
	{
		uint16_t mask = 0;
		assert(str[i] >= 0 && str[i] <= 3);
		mask = str[i];
		motif <<= 2;
		motif |= mask;
	}
	CharType code = getNumBefore(motifLength) + motif * motifLength + overhang;
	assert((end - start) / motifLength < std::numeric_limits<LengthType>::max()-1);
	LengthType runLength = ((end - start) - overhang) / motifLength;
	assert(runLength > 0);
	assert((((end - start) - overhang) % motifLength) == 0);
	return std::make_pair(code, runLength);
}

std::string getSequence(CharType code, size_t runLength)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	code -= getNumBefore(motifLength);
	uint16_t motif = code / motifLength;
	uint16_t overhang = code % motifLength;
	std::string motifStr;
	for (size_t i = 0; i < motifLength; i++)
	{
		uint16_t charCode = (motif >> ((motifLength - i - 1) * 2)) & 3;
		switch(charCode)
		{
			case 0:
				motifStr += 'A';
				break;
			case 1:
				motifStr += 'C';
				break;
			case 2:
				motifStr += 'G';
				break;
			case 3:
				motifStr += 'T';
				break;
		}
	}
	std::string result;
	for (size_t i = 0; i < runLength; i++)
	{
		result += motifStr;
	}
	if (overhang > 0) result += motifStr.substr(0, overhang);
	return result;
}

std::string getExpandedSequence(const SequenceCharType& codes, const SequenceLengthType& lens)
{
	assert(codes.size() == lens.size());
	std::string result;
	for (size_t i = 0; i < codes.size(); i++)
	{
		result += getSequence(codes[i], lens[i]);
	}
	return result;
}

std::pair<SequenceCharType, SequenceLengthType> multiRLECompressOne(const SequenceCharType& str, const size_t maxMaskLength)
{
	assert(maxMaskLength <= MaxMotifLength);
	std::vector<std::tuple<size_t, size_t, uint8_t>> runs;
	for (size_t i = 0; i < str.size(); i++)
	{
		size_t j = i;
		while (j < str.size() && str[j] == str[i]) j += 1;
		if (runs.size() == 0 || std::get<0>(runs.back()) > i || std::get<1>(runs.back()) < j)
		{
			runs.emplace_back(i, j, 1);
		}
		for (size_t motifLength = 2; motifLength <= maxMaskLength; motifLength++)
		{
			if (i + motifLength * 2 > str.size()) break;
			bool hasRepeat = true;
			for (size_t j = 0; j < motifLength; j++)
			{
				if (str[i+j] != str[i+j+motifLength])
				{
					hasRepeat = false;
					break;
				}
			}
			if (!hasRepeat) continue;
			size_t runLength = 1;
			while (i+motifLength*runLength+motifLength <= str.size() && SequenceCharType { str.begin() + i, str.begin() + i + motifLength } == SequenceCharType { str.begin() + i + motifLength * runLength, str.begin() + i + motifLength * runLength + motifLength })
			{
				runLength += 1;
			}
			assert(runLength >= 2);
			size_t overhang = 0;
			for (size_t j = 0; j < motifLength; j++)
			{
				if (i+motifLength*runLength+j >= str.size()) break;
				if (str[i+j] != str[i+motifLength*runLength+j]) break;
				overhang += 1;
			}
			assert(overhang < motifLength);
			size_t lengthHere = motifLength*runLength+overhang;
			assert(runs.size() > 0);
			if (std::get<0>(runs.back()) <= i && std::get<1>(runs.back()) >= i + lengthHere)
			{
				continue;
			}
			assert(runs.size() == 0 || (i > std::get<0>(runs.back()) || i + lengthHere > std::get<1>(runs.back())));
			while (runs.size() > 0 && std::get<0>(runs.back()) == i && std::get<1>(runs.back()) < i + lengthHere)
			{
				runs.pop_back();
			}
			assert(runs.size() == 0 || (i > std::get<0>(runs.back()) && i + lengthHere > std::get<1>(runs.back())));
			runs.emplace_back(i, i + lengthHere, motifLength);
		}
	}
	std::vector<std::tuple<size_t, size_t, uint8_t>> nonOverlappingRuns;
	size_t lastOneChar = 0;
	for (size_t i = 1; i < runs.size(); i++)
	{
		assert(std::get<0>(runs[i]) > std::get<0>(runs[i-1]));
		assert(std::get<1>(runs[i]) > std::get<1>(runs[i-1]));
		if (lastOneChar > std::get<0>(runs[i]))
		{
			for (size_t j = lastOneChar; j < std::get<1>(runs[i-1]); j++)
			{
				// maybe todo: hpc in overlapping parts
				nonOverlappingRuns.emplace_back(j, j+1, 1);
			}
			lastOneChar = std::get<1>(runs[i-1]);
			continue;
		}
		assert(lastOneChar <= std::get<0>(runs[i]));
		if (std::get<0>(runs[i]) >= std::get<1>(runs[i-1]))
		{
			assert(std::get<0>(runs[i]) == std::get<1>(runs[i-1]));
			assert(std::max(lastOneChar, std::get<0>(runs[i-1])) < std::get<1>(runs[i-1]));
			nonOverlappingRuns.emplace_back(std::max(lastOneChar, std::get<0>(runs[i-1])), std::get<1>(runs[i-1]), std::get<2>(runs[i-1]));
			continue;
		}
		assert(std::get<0>(runs[i]) < std::get<1>(runs[i-1]));
		assert(std::max(lastOneChar, std::get<0>(runs[i-1])) <= std::get<0>(runs[i]));
		if (std::max(lastOneChar, std::get<0>(runs[i-1])) < std::get<0>(runs[i])) nonOverlappingRuns.emplace_back(std::max(lastOneChar, std::get<0>(runs[i-1])), std::get<0>(runs[i]), std::get<2>(runs[i-1]));
		lastOneChar = std::get<1>(runs[i-1]);
		for (size_t j = std::get<0>(runs[i]); j < std::get<1>(runs[i-1]); j++)
		{
			// maybe todo: hpc in overlapping parts
			nonOverlappingRuns.emplace_back(j, j+1, 1);
		}
	}
	if (lastOneChar > std::get<0>(runs.back()))
	{
		if (lastOneChar != std::get<1>(runs.back())) nonOverlappingRuns.emplace_back(lastOneChar, std::get<1>(runs.back()), std::get<2>(runs.back()));
	}
	else
	{
		nonOverlappingRuns.emplace_back(runs.back());
	}
	assert(std::get<0>(nonOverlappingRuns[0]) == 0);
	assert(std::get<1>(nonOverlappingRuns.back()) == str.size());
	std::pair<SequenceCharType, SequenceLengthType> result;
	for (size_t i = 0; i < nonOverlappingRuns.size(); i++)
	{
		assert(i == 0 || std::get<1>(nonOverlappingRuns[i-1]) == std::get<0>(nonOverlappingRuns[i]));
		assert(std::get<1>(nonOverlappingRuns[i]) > std::get<0>(nonOverlappingRuns[i]));
		CharType code;
		LengthType runLength;
		uint8_t motifLength = std::get<2>(nonOverlappingRuns[i]);
		if (motifLength > std::get<1>(nonOverlappingRuns[i]) - std::get<0>(nonOverlappingRuns[i])) motifLength = std::get<1>(nonOverlappingRuns[i]) - std::get<0>(nonOverlappingRuns[i]);
		std::tie(code, runLength) = getCodeAndRunlength(str, std::get<0>(nonOverlappingRuns[i]), std::get<1>(nonOverlappingRuns[i]), motifLength);
		result.first.emplace_back(code);
		result.second.emplace_back(runLength);
	}
	return result;
}

std::vector<std::pair<SequenceCharType, SequenceLengthType>> multiRLECompress(const SequenceCharType& str, const size_t maxMaskLength)
{
	assert(maxMaskLength <= MaxMotifLength);
	std::vector<std::pair<SequenceCharType, SequenceLengthType>> result;
	size_t lastBreak = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		switch(str[i])
		{
			case 0:
			case 1:
			case 2:
			case 3:
				break;
			default:
				if (i > lastBreak+1) result.emplace_back(multiRLECompressOne(SequenceCharType { str.begin() + lastBreak, str.begin() + i }, maxMaskLength));
				lastBreak = i;
				break;
		}
	}
	if (lastBreak != str.size()-1) result.emplace_back(multiRLECompressOne(SequenceCharType { str.begin() + lastBreak, str.end() }, maxMaskLength));
	return result;
}

SequenceCharType revCompRLE(const SequenceCharType& codes)
{
	SequenceCharType result;
	result.resize(codes.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = complement(codes[codes.size()-1-i]);
	}
	return result;
}

size_t getExpandedLength(CharType code, LengthType length)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	code -= getNumBefore(motifLength);
	uint16_t overhang = code % motifLength;
	return length * motifLength + overhang;
}

size_t getOverlapFromRLE(const std::vector<std::pair<SequenceCharType, SequenceLengthType>>& unitigSequences, std::pair<size_t, bool> fromUnitig, size_t rleOverlap)
{
	if (fromUnitig.second)
	{
		SequenceCharType str { unitigSequences[fromUnitig.first].first.end() - rleOverlap, unitigSequences[fromUnitig.first].first.end() };
		SequenceLengthType lens { unitigSequences[fromUnitig.first].second.end() - rleOverlap, unitigSequences[fromUnitig.first].second.end() };
		return getExpandedSequence(str, lens).size();
	}
	else
	{
		SequenceCharType str { unitigSequences[fromUnitig.first].first.begin(), unitigSequences[fromUnitig.first].first.begin() + rleOverlap };
		str = revCompRLE(str);
		SequenceLengthType lens { unitigSequences[fromUnitig.first].second.rend() - rleOverlap, unitigSequences[fromUnitig.first].second.rend() };
		return getExpandedSequence(str, lens).size();
	}
	// assert(fromUnitig.first < unitigSequences.size());
	// assert(rleOverlap < unitigSequences[fromUnitig.first].first.size());
	// size_t result = 0;
	// for (size_t i = 0; i < rleOverlap; i++)
	// {
	// 	size_t index = i;
	// 	if (fromUnitig.second) index = unitigSequences[fromUnitig.first].first.size() - i - 1;
	// 	unsigned char chr = unitigSequences[fromUnitig.first].second[index];
	// 	if (chr >= 0 && chr <= 4)
	// 	{
	// 		result += unitigSequences[fromUnitig.first].second[index];
	// 	}
	// 	else if (chr >= 5 && chr <= 16)
	// 	{
	// 		result += unitigSequences[fromUnitig.first].second[index] * 2 + 1;
	// 	}
	// 	else if (chr >= 17 && chr <= 28)
	// 	{
	// 		result += unitigSequences[fromUnitig.first].second[index] * 2 ;
	// 	}
	// 	else
	// 	{
	// 		assert(false);
	// 	}
	// }
	// return result;
}

std::vector<size_t> getRLEExpandedPositions(const SequenceCharType& seq, const SequenceLengthType& lens)
{
	std::vector<size_t> result;
	result.resize(seq.size()+1);
	result[0] = 0;
	for (size_t i = 0; i < seq.size(); i++)
	{
		result[i+1] = result[i] + getExpandedLength(seq[i], lens[i]);
	}
	return result;
}
