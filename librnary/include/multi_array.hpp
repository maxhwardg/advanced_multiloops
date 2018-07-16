//
// Created by max on 6/23/16.
// Contains implementations of multi-dimensional arrays optimized for RNA folding DP tables.

#ifndef RNARK_MULTI_ARRAY_HPP
#define RNARK_MULTI_ARRAY_HPP

#include <vector>
#include <array>
#include <cstdarg>

namespace librnary {
/**
 * A multidimensional array optimized for the kinds of dynamic programming algorithms used in RNA research.
 * These typically have a large number of dimensions, each of which only has a small-moderate number of indexes.
 * This class, although very simple, can be a lot faster than nested vectors if pointer chasing is causing a problem.
 * Let's assume a random-like access pattern. Now, say you a 50*50*50 array, the MultiArray should be a lot faster.
 * The hardcoded ArrayX versions are faster, however. MultiArray should only be used for higher dimensional tables.
 */
template<typename T, unsigned NDIMS>
class MultiArray {
protected:
	std::array<unsigned, NDIMS> dims;
	std::vector<T> elems;
public:
	MultiArray() = default;
	MultiArray(const std::array<unsigned, NDIMS> &_dims, const T &initv)
		: dims(_dims) {
		unsigned sz = 1;
		for (unsigned d : dims) {
			sz *= d;
		}
		elems.assign(sz, initv);
	}
	template<typename... Args>
	T &operator()(Args... args) {
		return (*this)[{args...}];
	}
	T &operator[](const std::array<int, NDIMS> &loc) {
		int i = 0, mult = 1;
		for (int d = (int) loc.size() - 1; d >= 0; --d) {
			i += mult * loc[d];
			mult *= dims[d];
		}
		return elems[i];
	}
};

template<typename T>
class Array1D {
	T *arr;
	size_t c;
	bool owner = false;
public:
	Array1D()
		: c(0) {}
	Array1D(T *_arr, size_t _c)
		: arr(_arr), c(_c) {}
	Array1D(size_t _c, T base_val)
		: c(_c) {
		owner = true;
		arr = new T[c];
		std::fill(arr, arr + c, base_val);
	}
	void operator=(const Array1D<T> &base) {
		if (owner)
			delete[] arr;
		c = base.c;
		if (base.owner) {
			owner = true;
			arr = new T[c];
			std::copy(base.arr, base.arr + c, arr);
		} else {
			arr = base.arr;
		}
	}
	Array1D(const Array1D<T> &base) {
		(*this) = base;
	}
	~Array1D() {
		if (owner)
			delete[] arr;
	}
	const T &operator[](size_t loc) const {
		assert(loc < c);
		return arr[loc];
	}

	T &operator[](size_t loc) {
		assert(loc < c);
		return arr[loc];
	}
};

template<typename T>
class Array2D {
	T *arr;
	size_t r, c;
	bool owner = false;
public:
	Array2D()
		: r(0) {}
	Array2D(T *_arr, size_t _r, size_t _c)
		: arr(_arr), r(_r), c(_c) {}
	Array2D(size_t _r, size_t _c, T base_val)
		: r(_r), c(_c) {
		arr = new T[r * c];
		std::fill(arr, arr + r * c, base_val);
		owner = true;
	}
	void operator=(const Array2D<T> &base) {
		if (owner)
			delete[] arr;
		c = base.c;
		r = base.r;
		if (base.owner) {
			owner = true;
			arr = new T[r*c];
			std::copy(base.arr, base.arr + r*c, arr);
		} else {
			arr = base.arr;
		}
	}
	Array2D(const Array2D<T> &base) {
		(*this) = base;
	}
	~Array2D() {
		if (owner)
			delete[] arr;
	}
	const Array1D<T> operator[](size_t loc) const {
		assert(loc < r);
		return Array1D<T>(arr + loc * c, c);
	}
	Array1D<T> operator[](size_t loc) {
		assert(loc < r);
		return Array1D<T>(arr + loc * c, c);
	}
};

template<typename T>
class Array3D {
	T *arr;
	size_t rr, r, c;
	bool owner = false;
public:
	Array3D()
		: rr(0) {}
	Array3D(T *_arr, size_t _rr, size_t _r, size_t _c)
		: arr(_arr), rr(_rr), r(_r), c(_c) {}
	Array3D(size_t _rr, size_t _r, size_t _c, T base_val)
		: rr(_rr), r(_r), c(_c) {
		arr = new T[rr * r * c];
		std::fill(arr, arr + rr * r * c, base_val);
		owner = true;
	}
	Array3D(const Array3D<T> &base) {
		(*this) = base;
	}
	void operator=(const Array3D<T> &base) {
		if (owner)
			delete[] arr;
		c = base.c;
		r = base.r;
		rr = base.rr;
		if (base.owner) {
			owner = true;
			arr = new T[rr * r * c];
			std::copy(base.arr, base.arr + rr * r * c, arr);
		} else {
			arr = base.arr;
		}
	}
	~Array3D() {
		if (owner)
			delete[] arr;
	}
	const Array2D<T> operator[](size_t loc) const {
		assert(loc < rr);
		return Array2D<T>(arr + loc * r * c, r, c);
	}
	Array2D<T> operator[](size_t loc) {
		assert(loc < rr);
		return Array2D<T>(arr + loc * r * c, r, c);
	}
};

template<typename T>
class Array4D {
	T *arr;
	size_t rrr, rr, r, c;
	bool owner = false;
public:
	Array4D()
		: rrr(0) {}
	Array4D(T *_arr, size_t _rrr, size_t _rr, size_t _r, size_t _c)
		: arr(_arr), rrr(_rrr), rr(_rr), r(_r), c(_c) {}
	Array4D(size_t _rrr, size_t _rr, size_t _r, size_t _c, T base_val)
		: rrr(_rrr), rr(_rr), r(_r), c(_c) {
		arr = new T[rrr * rr * r * c];
		std::fill(arr, arr + rrr * rr * r * c, base_val);
		owner = true;
	}
	Array4D(const Array4D<T> &base) {
		(*this) = base;
	}
	void operator=(const Array4D<T> &base) {
		if (owner)
			delete[] arr;
		c = base.c;
		r = base.r;
		rr = base.rr;
		rrr = base.rrr;
		if (base.owner) {
			owner = true;
			arr = new T[rrr * rr * r * c];
			std::copy(base.arr, base.arr + rrr * rr * r * c, arr);
		} else {
			arr = base.arr;
		}
	}
	~Array4D() {
		if (owner)
			delete[] arr;
	}
	Array3D<T> operator[](size_t loc) {
		assert(loc < rrr);
		return Array3D<T>(arr + loc * rr * r * c, rr, r, c);
	}
	const Array3D<T> operator[](size_t loc) const {
		assert(loc < rrr);
		return Array3D<T>(arr + loc * rr * r * c, rr, r, c);
	}
};

}

#endif //RNARK_MULTI_ARRAY_HPP
