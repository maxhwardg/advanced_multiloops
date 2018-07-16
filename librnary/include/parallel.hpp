//
// Created by max on 7/27/17.
//

#ifndef RNARK_PARALLEL_HPP
#define RNARK_PARALLEL_HPP

#include <iterator>
#include <thread>
#include <future>
#include <vector>
#include <algorithm>
#include <cmath>

namespace librnary {

/**
 * Transforms one vector into another. Similar to (but not the same as) std::transform.
 * @tparam T1 Type of elements in the source vector.
 * @tparam T2 Type of elements in the target vector.
 * @tparam Func Transform function type.
 * @param A Source vector.
 * @param B Target vector.
 * @param f Tranformation function.
 * @param threads Number of threads to use.
 */
template<typename T1, typename T2, typename Func>
void parallel_transform(const std::vector<T1> &A,
						std::vector<T2> &B,
						const Func &f,
						size_t threads = std::thread::hardware_concurrency()) {
	using namespace std;
	auto block_sz = static_cast<size_t>(ceil(A.size() / static_cast<double>(threads)));
	auto execute = [&](size_t iA, size_t jA) {
		for (size_t i = iA; i < jA; ++i) {
			B[i] = f(A[i]);
		}
	};
	vector<thread> thread_list;
	for (size_t i = 0; i < threads; ++i) {
		thread_list.emplace_back(execute, i * block_sz, min(A.size(), (i + 1) * block_sz));
	}
	for (size_t i = 0; i < threads; ++i) {
		thread_list[i].join();
	}
}

namespace {
template<typename It>
It p_max_helper(It begin, It end, size_t threads) {
	if (threads == 1 || std::distance(begin, end) <= 1)
		return std::max_element(begin, end);
	else {
		It l, r;
		std::thread left([&]() {
			l = p_max_helper(begin, begin+std::distance(begin, end)/2, threads/2);
		});
		std::thread right([&]() {
			r = p_max_helper(begin+std::distance(begin, end)/2, end, threads/2);
		});
		left.join();
		right.join();
		if (*r > *l)
			return r;
		return l;
	}
}
template<typename It>
It p_min_helper(It begin, It end, size_t threads) {
	if (threads == 1 || std::distance(begin, end) <= 1)
		return std::min_element(begin, end);
	else {
		It l, r;
		std::thread left([&]() {
			l = p_min_helper(begin, begin+std::distance(begin, end)/2, threads/2);
		});
		std::thread right([&]() {
			r = p_min_helper(begin+std::distance(begin, end)/2, end, threads/2);
		});
		left.join();
		right.join();
		if (*r < *l)
			return r;
		return l;
	}
}
}

/**
 * The same as to std::max_element.
 * @tparam It Iterator type.
 * @param begin Begin iterator.
 * @param end End iterator.
 * @return Iterator to the maximum element.
 */
template<typename It>
It parallel_max_element(It begin, It end, size_t threads = std::thread::hardware_concurrency()) {
	assert(threads >= 1);
	return p_max_helper(begin, end, threads);
}

/**
 * The same as to std::min_element.
 * @tparam It Iterator type.
 * @param begin Begin iterator.
 * @param end End iterator.
 * @return Iterator to the maximum element.
 */
template<typename It>
It parallel_min_element(It begin, It end, size_t threads = std::thread::hardware_concurrency()) {
	assert(threads >= 1);
	return p_min_helper(begin, end, threads);
}
}

#endif //RNARK_PARALLEL_HPP
