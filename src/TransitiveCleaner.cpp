#include "TransitiveCleaner.h"
#include "FastHasher.h"

TransitiveCleaner::TransitiveCleaner(size_t kmerSize, const HashList& hashlist)
{
	transitiveMiddle.resize(hashlist.size());
	getMiddles(kmerSize, hashlist);
}

std::vector<std::pair<size_t, bool>> TransitiveCleaner::insertMiddles(std::vector<std::pair<size_t, bool>> raw) const
{
	std::vector<std::pair<size_t, bool>> result;
	while (raw.size() >= 2)
	{
		auto from = raw[raw.size()-2];
		auto to = raw[raw.size()-1];
		if (transitiveMiddle[from].count(to) == 0)
		{
			result.push_back(raw.back());
			raw.pop_back();
			continue;
		}
		raw.pop_back();
		assert(transitiveMiddle[from].at(to).size() > 0);
		raw.insert(raw.end(), transitiveMiddle[from].at(to).begin(), transitiveMiddle[from].at(to).end());
		raw.push_back(to);
	}
	result.push_back(raw.back());
	raw.pop_back();
	std::reverse(result.begin(), result.end());
	return result;
}

void TransitiveCleaner::addMiddles(size_t kmerSize, std::pair<size_t, bool> start, std::pair<size_t, bool> end, LazyString& seq, const HashList& list, const std::unordered_set<size_t>& minimizerPrefixes)
{
	std::vector<std::pair<size_t, bool>> path;
	std::pair<size_t, bool> old = start;
	size_t oldpos = 0;

	uint64_t fwHash = start.second ? list.fakeFwHashes[start.first] : list.fakeBwHashes[start.first];
	uint64_t bwHash = start.second ? list.fakeBwHashes[start.first] : list.fakeFwHashes[start.first];
	FastHasher fwkmerHasher { kmerSize, fwHash, bwHash };
	assert(minimizerPrefixes.count(fwkmerHasher.hash()) == 1);
	for (size_t i = 1; i < seq.size() - kmerSize; i++)
	{
		fwkmerHasher.addChar(seq[i + kmerSize - 1]);
		fwkmerHasher.removeChar(seq[i - 1]);
		size_t hash = fwkmerHasher.hash();
		if (minimizerPrefixes.count(hash) == 1)
		{
			auto here = list.getNodeOrNull(seq.view(i, kmerSize));
			if (here.first == std::numeric_limits<size_t>::max())
			{
				continue;
			}
			path.push_back(here);
			std::pair<size_t, bool> canonFrom;
			std::pair<size_t, bool> canonTo;
			std::tie(canonFrom, canonTo) = canon(old, here);
			newSequenceOverlaps.emplace_back(canonFrom, canonTo, kmerSize - (i - oldpos));
			old = here;
			oldpos = i;
		}
	}

	if (path.size() > 0)
	{
		assert(old != start);
		std::pair<size_t, bool> canonFrom;
		std::pair<size_t, bool> canonTo;
		std::tie(canonFrom, canonTo) = canon(old, end);
		newSequenceOverlaps.emplace_back(canonFrom, canonTo, kmerSize - (seq.size() - kmerSize - oldpos));
		transitiveMiddle[start][end] = path;
	}
	else
	{
		assert(old == start);
	}
}

std::unordered_set<size_t> TransitiveCleaner::getMinimizerPrefixes(size_t kmerSize, const HashList& hashlist)
{
	std::unordered_set<size_t> result;
	result.insert(hashlist.fakeFwHashes.begin(), hashlist.fakeFwHashes.end());
	result.insert(hashlist.fakeBwHashes.begin(), hashlist.fakeBwHashes.end());
	return result;
}

void TransitiveCleaner::getMiddles(size_t kmerSize, const HashList& hashlist)
{
	std::unordered_set<size_t> minimizerPrefixes = getMinimizerPrefixes(kmerSize, hashlist);
	for (size_t i = 0; i < hashlist.sequenceOverlap.size(); i++)
	{
		std::pair<size_t, bool> fw { i, true };
		if (hashlist.sequenceOverlap[fw].size() > 0)
		{
			TwobitView seq = hashlist.getHashSequenceRLE(i);
			for (auto pair : hashlist.sequenceOverlap[fw])
			{
				assert(pair.first.first >= i);
				TwobitView second { pair.first.second ? hashlist.getHashSequenceRLE(pair.first.first) : hashlist.getRevCompHashSequenceRLE(pair.first.first) };
				LazyString lazy { seq, second, pair.second };
				addMiddles(kmerSize, fw, pair.first, lazy, hashlist, minimizerPrefixes);
			}
		}
		std::pair<size_t, bool> bw { i, false };
		if (hashlist.sequenceOverlap[bw].size() > 0)
		{
			TwobitView seq = hashlist.getRevCompHashSequenceRLE(i);
			for (auto pair : hashlist.sequenceOverlap[bw])
			{
				assert(pair.first.first >= i);
				TwobitView second { pair.first.second ? hashlist.getHashSequenceRLE(pair.first.first) : hashlist.getRevCompHashSequenceRLE(pair.first.first) };
				LazyString lazy { seq, second, pair.second };
				addMiddles(kmerSize, bw, pair.first, lazy, hashlist, minimizerPrefixes);
			}
		}
	}
}

