#ifndef GEOM_CONST_H
#define GEOM_CONST_H

#include <cstdint>

/*
    This file defines geometric constants related to geometry.
    It includes:
    - GRID: A scaling factor to convert between double precision coordinates and int64_t for tolerance handling.
    - TOLERANCE, LARGE_TOLERANCE, SAME_POINT: Constants defining various levels of precision for geometric calculations.
    - EPS: A small epsilon value for floating-point comparisons.
    - MIDDLE_POINT_SHIFT: A small shift value used when calculating middle points to avoid numerical issues.
    - EPS_INT64: An integer representation of the epsilon value for use in int64_t calculations.
*/

#ifndef M_PI
    constexpr double M_PI = 3.14159265358979323846;
#endif

namespace geom {
    constexpr int64_t GRID              = 100000000;    // This is the scale to be able to preserve the numerical tolerance
    constexpr int64_t TOLERANCE         = 1000;         // 1E-5 precision
    constexpr int64_t LARGE_TOLERANCE   = 10000;        // 1E-3 precision
    constexpr int64_t SAME_POINT        = 100;          // 1E-6 precision
    constexpr double EPS = 1e-8;
    constexpr double MIDDLE_POINT_SHIFT = 1e-6;
    constexpr int64_t EPS_INT64 = 1000;
}

#endif // GEOM_CONST_H