#include "pch.h"
#include "ParticleGrid.h"
#include "ParticleSystem.h"
#include <limits>

using Fluid::Particle;
using Fluid::ParticleGrid;
using Fluid::ParticleSystem;
ParticleGrid::ParticleGrid(const ParticleSystem* system):
    m_system(system)
{
}


ParticleGrid::~ParticleGrid()
{
}

void ParticleGrid::UpdateGrid()
{
    m_grid.clear();
    DirectX::XMVECTOR boxmin =
        DirectX::XMVectorReplicate(std::numeric_limits<float>::max());
    DirectX::XMVECTOR boxmax =
        DirectX::XMVectorReplicate(-std::numeric_limits<float>::max());
    for (auto& p : m_system->GetParticles())
    {
        DirectX::XMVECTOR position = XMLoadFloat3(&p.position);
        boxmin = DirectX::XMVectorMin(boxmin, position);
        boxmax = DirectX::XMVectorMax(boxmax, position);
    }
    DirectX::XMVECTOR boxSize =
        DirectX::XMVectorSubtract(boxmax, boxmin);
    DirectX::XMVECTOR cellNum =
        DirectX::XMVectorScale(boxSize, 1 / m_radius);
    cellNum = DirectX::XMVectorCeiling(cellNum);
    DirectX::XMStoreSInt3(&m_cellNum, cellNum);
    DirectX::XMStoreFloat3(&m_boxmin, boxmin);
    m_grid.resize((m_cellNum.x + 1) * (m_cellNum.y + 1) * (m_cellNum.z + 1), std::vector<int>());
    for (int i = 0; i < static_cast<int>(m_system->GetParticles().size()); ++i)
    {
        DirectX::XMVECTOR position =
            DirectX::XMLoadFloat3(&m_system->GetParticles()[i].position);
        DirectX::XMINT3 pos;
        DirectX::XMStoreSInt3(&pos, GetCellIndex(position));
        GetCell(pos.x, pos.y, pos.z).push_back(i);
    }
}

std::vector<int> ParticleGrid::GetNeighbors(int index)
{
    std::vector<int> result;
    DirectX::XMVECTOR point = DirectX::XMLoadFloat3(&m_system->GetParticles()[index].position);
    DirectX::XMINT3 cellIndex;
    DirectX::XMStoreSInt3(&cellIndex,
        GetCellIndex(point));
    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <= 1; ++k)
            {
                int x = cellIndex.x + i;
                int y = cellIndex.y + j;
                int z = cellIndex.z + k;
                if (x < 0 || x >= m_cellNum.x) continue;
                if (y < 0 || y >= m_cellNum.y) continue;
                if (z < 0 || z >= m_cellNum.z) continue;
                for (int pointId : GetCell(cellIndex.x + i, cellIndex.y + j, cellIndex.z + k))
                {
                    if (pointId == index) continue;
                    DirectX::XMVECTOR point2 =
                        DirectX::XMLoadFloat3(&m_system->GetParticles()[pointId].position);
                    DirectX::XMVECTOR distanceVector =
                        DirectX::XMVector3Length(
                            DirectX::XMVectorSubtract(point2, point));
                    float distance;
                    DirectX::XMStoreFloat(&distance, distanceVector);
                    if (distance < m_radius)
                    {
                        result.push_back(pointId);
                    }
                }
            }
    return result;
}
