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
private:
	std::vector<SmallType> smallIncreases;
	phmap::flat_hash_map<size_t, size_t> bigIncreases;
	mutable size_t lastIndex;
	mutable size_t lastValue;
};

#endif
