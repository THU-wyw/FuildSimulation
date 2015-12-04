#pragma once
#ifndef SPH_KERNEL_FUNCTIONS
#define SPH_KERNEL_FUNCTIONS
#include <amp.h>
#include <amp_math.h>
#include <amp_graphics.h>
#include <DirectXMath.h>
#include "AMPMathExtend.h"
namespace Fluid
{
    class Poly6Kernel
    {
    public:
        Poly6Kernel(float kernel) : KERNEL(kernel) {}
        inline float Calculate(Concurrency::graphics::float_3 ri, Concurrency::graphics::float_3 rj) const restrict(amp, cpu)
        {
            float r2 = Float3LengthSq(ri - rj);
            return VALUE * (KERNEL2 - r2) * (KERNEL2 - r2) * (KERNEL2 - r2);
        }

        inline Concurrency::graphics::float_3 CalculateGrad(
            Concurrency::graphics::float_3 ri,
            Concurrency::graphics::float_3 rj) const restrict(amp)
        {
            Concurrency::graphics::float_3 rij = ri - rj;
            float r2 = Float3LengthSq(rij);
            return GRAD_VALUE * (KERNEL2 - r2) * (KERNEL2 - r2)  * rij;
        }

        inline float CalculateLplc(
            Concurrency::graphics::float_3 ri,
            Concurrency::graphics::float_3 rj) const restrict(amp)
        {
            Concurrency::graphics::float_3 rij = ri - rj;
            float r2 = Float3LengthSq(rij);
            return LPLC_VALUE * (KERNEL2 - r2) * (r2 - 0.75f * (KERNEL2 - r2));
        }
    private:
        const float KERNEL;
        const float KERNEL2 = KERNEL * KERNEL;
        const float VALUE = 315.0f / (64.0f * DirectX::XM_PI * pow(KERNEL, 9));
        const float GRAD_VALUE = -945.0f / (32.0f * DirectX::XM_PI * pow(KERNEL, 9));
        const float LPLC_VALUE = 945.0f / (8.0f * DirectX::XM_PI * pow(KERNEL, 9));
    };

    class SpikyKernel
    {
    public:
        SpikyKernel(float kernel) : KERNEL(kernel) {}
        inline Concurrency::graphics::float_3 CalculateGrad(
            Concurrency::graphics::float_3 ri,
            Concurrency::graphics::float_3 rj) const restrict(amp)
        {
            Concurrency::graphics::float_3 rij = ri - rj;
            float r = Float3Length(rij);
            return r > 0.001f
                ? rij * GRAD_VALUE * (KERNEL - r) * (KERNEL - r) / r
                : Concurrency::graphics::float_3(0);
        }
    private:
        const float KERNEL;
        const float KERNEL2 = KERNEL * KERNEL;
        
        const float GRAD_VALUE = -45.0f / (DirectX::XM_PI * pow(KERNEL, 6));
    };

    class ViscosityKernel
    {
    public:
        ViscosityKernel(float kernel) : KERNEL(kernel) {}
        inline float CalculateLplc(Concurrency::graphics::float_3 ri, Concurrency::graphics::float_3 rj) const restrict(amp, cpu)
        {
            float r = Float3Length(ri - rj);
            return LPLC_VALUE * (KERNEL - r);
        }
    private:
        const float KERNEL;
        const float KERNEL2 = KERNEL * KERNEL;
        const float LPLC_VALUE = 45.0f / (DirectX::XM_PI * pow(KERNEL, 6));
    };
}

#endif