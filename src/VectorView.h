#ifndef VectorView_h
#define VectorView_h

#include <string_view>
#include <vector>

template <typename T>
class VectorView
{
public:
	VectorView(std::vector<T>& data, size_t startpos, size_t endpos) :
	data(data),
	startpos(startpos),
	endpos(endpos)
	{
	}
	size_t size() const
	{
		return endpos-startpos;
	}
	const T& operator[](size_t pos) const
	{
		return data[startpos+pos];
	}
	T& operator[](size_t pos)
	{
		return data[startpos+pos];
	}
	T* begin()
	{
		return data.data() + startpos;
	}
	T* end()
	{
		return data.data() + endpos;
	}
	std::string_view getView()
	{
		return std::string_view { (char*)(data.data() + startpos), size() * sizeof(T) };
	}
private:
	std::vector<T>& data;
	size_t startpos;
	size_t endpos;
};

#endif
