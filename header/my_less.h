#ifndef MYLESS_HPP
#define MYLESS_HPP
#include <queue>
#include <functional>
#include <type_traits>
#include "my_greater.h"

#if __GNUC__ >= 3
#include <ext/functional>
using __gnu_cxx::select2nd;
using __gnu_cxx::identity;
#else
#include "my_select.h"
#endif

template<class T, class S = my_identity<T> >
class my_less {
private:
	S _s;
public:
	inline bool operator()(T const &a, T const &b) const {
		return _s(a) < _s(b);
	}
};


#endif
