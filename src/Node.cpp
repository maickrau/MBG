#include <cassert>
#include "Node.h"

Node::Node()
{
	val = std::numeric_limits<size_t>::max();
}

Node::Node(size_t id, bool forward)
{
	assert(id < std::numeric_limits<size_t>::max() / 2);
	val = (id << 1) + (forward ? 1 : 0);
}

Node::Node(std::pair<size_t, bool> p)
{
	assert(p.first < std::numeric_limits<size_t>::max() / 2);
	val = (p.first << 1) + (p.second ? 1 : 0);
}

size_t Node::id() const
{
	assert(val != std::numeric_limits<size_t>::max());
	return val >> 1;
}

bool Node::forward() const
{
	assert(val != std::numeric_limits<size_t>::max());
	return (val & 1) == 1;
}

Node::operator std::pair<size_t, bool>() const
{
	assert(val != std::numeric_limits<size_t>::max());
	return std::make_pair(id(), forward());
}

bool Node::operator==(const std::pair<size_t, bool>& p) const
{
	return (id() == p.first) && (forward() == p.second);
}

bool Node::operator<(const Node& other) const
{
	return val < other.val;
}
