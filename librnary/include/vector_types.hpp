//
// Created by max on 7/4/16.
// Vector type shortcut. Useful for declaring DP arrays.

#ifndef RNARK_VECTOR_TYPES_HPP
#define RNARK_VECTOR_TYPES_HPP

#include <vector>

#include "energy.hpp"

namespace librnary {
template<typename T>
using V = std::vector<T>;
template<typename T>
using VV = std::vector<V<T>>;
template<typename T>
using VVV = std::vector<VV<T>>;
template<typename T>
using _4DV = std::vector<VVV<T>>;
template<typename T>
using _5DV = std::vector<_4DV<T>>;
template<typename T>
using _6DV = std::vector<_5DV<T>>;
template<typename T>
using _7DV = std::vector<_6DV<T>>;

typedef V<int> VI;
typedef VV<int> VVI;
typedef VVV<int> VVVI;
typedef _4DV<int> _4DVI;
typedef _5DV<int> _5DVI;
typedef _6DV<int> _6DVI;
typedef _7DV<int> _7DVI;

typedef V<energy_t> VE;
typedef VV<energy_t> VVE;
typedef VVV<energy_t> VVVE;
typedef _4DV<energy_t> _4DVE;
typedef _5DV<energy_t> _5DVE;
typedef _6DV<energy_t> _6DVE;
typedef _7DV<energy_t> _7DVE;
}

#endif //RNARK_VECTOR_TYPES_HPP
