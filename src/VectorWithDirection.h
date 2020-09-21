#ifndef VectorWithDirection_h
#define VectorWithDirection_h

template <typename T>
class VectorWithDirection
{
public:
	void resize(size_t newSize)
	{
		forward.resize(newSize);
		backward.resize(newSize);
	}
	void resize(size_t newSize, T item)
	{
		forward.resize(newSize, item);
		backward.resize(newSize, item);
	}
	size_t size() const
	{
		return forward.size();
	}
	template <typename... Args>
	void emplace_back(Args... args)
	{
		forward.emplace_back(args...);
		backward.emplace_back(args...);
	}
	typename std::vector<T>::reference at(std::pair<size_t, bool> index)
	{
		assert(index.first < forward.size());
		return (*this)[index];
	}
	typename std::vector<T>::const_reference at(std::pair<size_t, bool> index) const
	{
		assert(index.first < forward.size());
		return (*this)[index];
	}
	typename std::vector<T>::reference operator[](std::pair<size_t, bool> index)
	{
		if (index.second) return forward[index.first];
		return backward[index.first];
	}
	typename std::vector<T>::const_reference operator[](std::pair<size_t, bool> index) const
	{
		if (index.second) return forward[index.first];
		return backward[index.first];
	}
	void clear()
	{
		forward.clear();
		backward.clear();
	}
private:
	std::vector<T> forward;
	std::vector<T> backward;
};

#endif
