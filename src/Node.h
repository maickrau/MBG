#ifndef Node_h
#define Node_h

#include <cstdlib>
#include <tuple>

namespace MBG
{

class Node
{
public:
	Node();
	Node(size_t id, bool forward);
	Node(std::pair<size_t, bool> p);
	size_t id() const;
	bool forward() const;
	operator std::pair<size_t, bool>() const;
	bool operator==(const std::pair<size_t, bool>& p) const;
	bool operator<(const Node& other) const;
private:
	size_t val;
};

}

#endif
