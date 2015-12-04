#pragma once
#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H
#include <DirectXMath.h>

namespace Fluid
{
    struct Particle
    {
        DirectX::XMFLOAT3 position;
        DirectX::XMFLOAT3 velocity;
        DirectX::XMFLOAT3 acceleration;
        DirectX::XMFLOAT3 velocityEval;

        float density;
        float pressure;
    };
}

#endif