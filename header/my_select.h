#ifndef my_select_hpp
#define my_select_hpp

template<typename Pair>
struct select2nd {

	const typename Pair::second_type& operator()(const Pair& x) const {
		return x.second;
	}
	typename Pair::second_type& operator()(Pair& x) {
		return x.second;
	}
};

//struct select2nd
//{
//	template< typename K, typename V >
//	const V& operator()(std::pair<K, V> const& p) const
//	{
//		return p.second;
//	}
//};

#endif