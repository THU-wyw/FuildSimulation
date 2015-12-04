#pragma once
#ifndef SPH_AMP_MATH_EXTEND
#define SPH_AMP_MATH_EXTEND

#include <amp_graphics.h>
#include <amp_math.h>

namespace Fluid
{
    template<typename T>
    constexpr inline T min(const T& a, const T& b) restrict(amp, cpu)
    {
        return a < b ? a : b;
    }

    template<typename T>
    constexpr inline T max(const T& a, const T& b) restrict(amp, cpu)
    {
        return a > b ? a : b;
    }

    inline float Float3Dot(
        const Concurrency::graphics::float_3& a,
        const Concurrency::graphics::float_3& b)
        restrict(amp, cpu)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline Concurrency::graphics::float_3 Float3Cross(
        const Concurrency::graphics::float_3& a,
        const Concurrency::graphics::float_3& b)
        restrict(amp, cpu)
    {
        return Concurrency::graphics::float_3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
    }

    inline float Float3Length(const Concurrency::graphics::float_3& x) restrict(amp, cpu)
    {
        return Concurrency::fast_math::sqrt(Float3Dot(x, x));
    }

    inline float Float3LengthSq(const Concurrency::graphics::float_3& x) restrict(amp, cpu)
    {
        return Float3Dot(x, x);
    }

    inline Concurrency::graphics::float_3 Float3Min(
        const Concurrency::graphics::float_3& a,
        const Concurrency::graphics::float_3& b)
        restrict(amp)
    {
        return Concurrency::graphics::float_3(
            Concurrency::fast_math::fminf(a.x, b.x),
            Concurrency::fast_math::fminf(a.y, b.y),
            Concurrency::fast_math::fminf(a.z, b.z)
            );
    }

    inline Concurrency::graphics::float_3 Float3Max(
        const Concurrency::graphics::float_3& a,
        const Concurrency::graphics::float_3& b)
        restrict(amp)
    {
        return Concurrency::graphics::float_3(
            Concurrency::fast_math::fmaxf(a.x, b.x),
            Concurrency::fast_math::fmaxf(a.y, b.y),
            Concurrency::fast_math::fmaxf(a.z, b.z)
            );
    }
}
#endif