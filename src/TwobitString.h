#ifndef TwobitString_h
#define TwobitString_h

#include <string>
#include <string_view>
#include <vector>

class TwobitString
{
public:
	TwobitString();
	TwobitString(const std::string_view& str);
	char get(size_t i) const;
	void set(size_t i, char c);
	size_t size() const;
	TwobitString& operator=(const std::string& str);
	TwobitString& operator=(const std::string_view& str);
	void resize(size_t size);
	void push_back(char c);
	template <typename Iter>
	void insert(Iter start, Iter end)
	{
		Iter i = start;
		while (i != end)
		{
			push_back(*i);
			++i;
		}
	}
private:
	std::vector<unsigned char> data;
	size_t realSize;
};


class TwobitView
{
public:
	TwobitView(const TwobitString& str, const size_t start, const size_t end);
	char operator[](size_t i) const;
	size_t size() const;
	std::string toString() const;
	std::string toSubstring(size_t substrStart) const;
private:
	const TwobitString& str;
	const size_t start;
	const size_t end;
};

#endif
