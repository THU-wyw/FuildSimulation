#pragma once
#ifndef SPH_PARTICLE_GRID_H
#define SPH_PARTICLE_GRID_H
#include <DirectXMath.h>
#include <DirectXCollision.h>
#include <vector>
#include "Particle.h"

namespace Fluid
{
    class ParticleSystem;
    class ParticleGrid
    {
    public:
        ParticleGrid(const ParticleSystem* system);
        ParticleGrid(const ParticleGrid&) = delete;
        ~ParticleGrid();
        void UpdateGrid();
        std::vector<int> GetNeighbors(int index);
        void SetRadius(float radius) { m_radius = radius; }
    private:
        const ParticleSystem* m_system;
        std::vector<std::vector<int>> m_grid;
        float m_radius;
        DirectX::XMINT3 m_cellNum;
        DirectX::XMFLOAT3 m_boxmin;

        inline std::vector<int>& GetCell(int x, int y, int z)
        {
            return m_grid[((x * m_cellNum.y) + y) * m_cellNum.z + z];
        }
        inline DirectX::XMVECTOR GetCellIndex(DirectX::XMVECTOR point)
        {
            DirectX::XMVECTOR offset =
                DirectX::XMVectorSubtract(
                    point,
                    DirectX::XMLoadFloat3(&m_boxmin));
            DirectX::XMVECTOR cellIndexVector =
                DirectX::XMVectorScale(offset, 1 / m_radius);
            return DirectX::XMVectorFloor(cellIndexVector);
        }
    };
}

#endif