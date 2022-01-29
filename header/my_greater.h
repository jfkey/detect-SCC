#ifndef MYGREATER_HPP
#define MYGREATER_HPP
#include <queue>
#include <functional>
#include <type_traits>

#if __GNUC__ >= 3
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::identity;
#else
#include "my_select.h"
#endif

template<class T, class S = my_identity<T> >
class my_greater {
private:
	S _s;
public:
	inline bool operator()(T const &a, T const &b) const {
		return _s(b) < _s(a);
	}
};



template <typename T>
struct my_identity
{
	T operator()(T x) const { return x; }
};

#endif
