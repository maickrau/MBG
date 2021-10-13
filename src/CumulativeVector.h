#ifndef CumulativeVector_h
#define CumulativeVector_h

template <typename SmallType>
class CumulativeVector
{
public:
	CumulativeVector() :
	smallIncreases(),
	bigIncreases(),
	lastIndex(0),
	lastValue(0)
	{}
	size_t size() const
	{
		return smallIncreases.size();
	}
	void push_back(size_t value)
	{
		if (size() == 0) lastValue = value;
		size_t oldValue = 0;
		if (size() > 0) oldValue = get(size() - 1);
		assert(value >= oldValue);
		if ((size_t)value - (size_t)oldValue < (size_t)std::numeric_limits<SmallType>::max())
		{
			smallIncreases.push_back((SmallType)((size_t)value - (size_t)oldValue));
		}
		else
		{
			bigIncreases[size()] = (size_t)value - (size_t)oldValue;
			smallIncreases.push_back(0);
		}
	}
	size_t operator[](size_t index) const
	{
		return get(index);
	}
	size_t back() const
	{
		assert(size() > 0);
		return get(size()-1);
	}
	void clear()
	{
		smallIncreases.clear();
		bigIncreases.clear();
		lastIndex = 0;
		lastValue = 0;
	}
	size_t get(size_t index) const
	{
		assert(index < size());
		if (index < lastIndex)
		{
			lastIndex = 0;
			lastValue = smallIncreases[0];
			auto found = bigIncreases.find(0);
			if (found != bigIncreases.end()) lastValue += found->second;
		}
		for (size_t i = lastIndex+1; i <= index; i++)
		{
			lastValue += smallIncreases[i];
			auto found = bigIncreases.find(i);
			if (found != bigIncreases.end()) lastValue += found->second;
		}
		lastIndex = index;
		return lastValue;
	}
	void eraseRange(size_t start, size_t end)
	{
		assert(start == 0 || end == size());
		if (start == 0 && end == size())
		{
			clear();
		}
		else if (start == 0)
		{
			assert(end < size());
			size_t newFirst = get(end);
			smallIncreases.erase(smallIncreases.begin(), smallIncreases.begin() + end);
			phmap::flat_hash_map<size_t, size_t> newBigs;
			for (auto pair : bigIncreases)
			{
				if (pair.first < end) continue;
				newBigs[pair.first - end] = pair.second;
			}
			std::swap(bigIncreases, newBigs);
			if (newFirst < (size_t)std::numeric_limits<SmallType>::max())
			{
				smallIncreases[0] = newFirst;
				if (bigIncreases.count(0) == 1) bigIncreases.erase(0);
			}
			else
			{
				smallIncreases[0] = 0;
				bigIncreases[0] = newFirst;
			}
		}
		else
		{
			assert(start > 0);
			smallIncreases.erase(smallIncreases.begin() + start, smallIncreases.end());
			phmap::flat_hash_map<size_t, size_t> newBigs;
			for (auto pair : bigIncreases)
			{
				if (pair.first >= start) continue;
				newBigs[pair.first] = pair.second;
			}
			std::swap(bigIncreases, newBigs);
		}
		lastIndex = 0;
		lastValue = get(0);
	}
	template <typename Iterator>
	void insertEnd(Iterator start, Iterator end)
	{
		while (start != end)
		{
			push_back(*start);
			++start;
		}
	}
private:
	std::vector<SmallType> smallIncreases;
	phmap::flat_hash_map<size_t, size_t> bigIncreases;
	mutable size_t lastIndex;
	mutable size_t lastValue;
};

#endif
