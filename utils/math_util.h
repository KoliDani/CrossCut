#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

/*
    Simple math utilities to avoid conflicts with Windows.h min/max macros.
    This file should be included after all Windows headers and before any code that uses min/max.
*/
namespace MathUtil {
    // Template min/max in global namespace for convenience
    template<typename T>
    inline const T& min(const T& a, const T& b) { return (b < a) ? b : a; }

    template<typename T>
    inline const T& max(const T& a, const T& b) { return (a < b) ? b : a; }
}

#endif // MATH_UTIL_H
