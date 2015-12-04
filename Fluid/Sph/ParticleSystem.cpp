#include "pch.h"
#include "ParticleSystem.h"
#include <limits>
#include <fstream>
#include "AMPMathExtend.h"
#include "AMPAlgorithmExtend.h"
#include "KernelFuntions.h"
#include "../Lib/amp_algorithms.h"
using Fluid::ParticleSystem;
using Concurrency::graphics::float_3;
using Concurrency::array;
using Concurrency::index;
using Concurrency::graphics::int_3;

ParticleSystem::ParticleSystem(ID3D11Device2* device) :
    m_particles(),
    m_grid(this),
    m_accelerator_view(Concurrency::direct3d::create_accelerator_view(device))
{
    m_grid.SetRadius(m_kernel);
}

ParticleSystem::~ParticleSystem()
{
}

void ParticleSystem::InitParticles()
{
    m_box_min = { 0, 0, 0 };
    m_box_max = { 0.64f, 0.64f, 0.64f };
    float spacing = 0.5f;
    for (float x = 0.0f; x < 0.64f * 0.6f; x += m_kernel / 2)
        for (float y = 0.0f; y < 0.64f * 0.9f; y += m_kernel / 2)
            for (float z = 0.0f; z < 0.64f * 0.6f; z += m_kernel / 2)
            {
                m_points_vector.push_back({ x, y, z });
                m_particles.push_back({
                    {x, y, z},
                    {0, 0, 0},
                    {0, 0, 0},
                    {0, 0, 0},
                    0,
                    0
                });
            }
    m_points = std::make_unique<array<float_3>>(m_points_vector.size(), m_points_vector.begin(), m_points_vector.end(), m_accelerator_view);
    m_velocity = std::make_unique<array<float_3>>(m_points->extent, m_accelerator_view);
    m_velocity_eval = std::make_unique<array<float_3>>(m_points->extent, m_accelerator_view);
    m_density = std::make_unique<array<float>>(m_points->extent, m_accelerator_view);
    m_pressure = std::make_unique<array<float>>(m_points->extent, m_accelerator_view);
    m_force = std::make_unique<array<float_3>>(m_points->extent, m_accelerator_view);
    m_points_next = std::make_unique<array<int>>(m_points->extent, m_accelerator_view);
    m_neighbor_indice = std::make_unique<array<int>>(m_points_vector.size() + 1, m_accelerator_view, Concurrency::access_type_read_write);
}

void ParticleSystem::Update()
{
    GenerateNeighborTable();
    //TestNeighbors();
    //std::vector<float_3> d = *m_points;
    //std::vector<int> a = *m_neighbors;
    ComputeDensity();
    //TestDensity();
    //std::vector<float> b = *m_density;
    //std::vector<float_3> f = *m_force;
    ComputeInternalForce();
    //TestForce();
    //std::vector<float_3> e = *m_force;
    Advect();
    //std::vector<float_3> c = *m_points;
}

void ParticleSystem::GetBox()
{
    m_cell_box_min = float_3(std::numeric_limits<float>::max());
    m_cell_box_max = float_3(std::numeric_limits<float>::lowest());
    m_points_vector = *m_points;
    for (auto p : m_points_vector)
    {
        m_cell_box_min.x = std::min(m_cell_box_min.x, p.x);
        m_cell_box_min.y = std::min(m_cell_box_min.y, p.y);
        m_cell_box_min.z = std::min(m_cell_box_min.z, p.z);
        m_cell_box_max.x = std::max(m_cell_box_max.x, p.x);
        m_cell_box_max.y = std::max(m_cell_box_max.y, p.y);
        m_cell_box_max.z = std::max(m_cell_box_max.z, p.z);
    }
    Concurrency::array_view<float_3> points_view(m_points_vector.size(), m_points_vector);
    float_3 temp_min = amp_algorithms::reduce(points_view,
        [](const float_3& a, const float_3& b) restrict(cpu, amp) -> float_3
    {
        return{ min(a.x, b.x), min(a.y, b.y), min(a.z, b.z) };
    });

    //temp_min = amp_stl_algorithms::reduce(amp_stl_algorithms::begin(points_view), 
    //    amp_stl_algorithms::end(points_view), float_3(std::numeric_limits<float>::max()),
    //    [](const float_3& a, const float_3& b) restrict(cpu, amp) -> float_3
    //{
    //    return{ min(a.x, b.x), min(a.y, b.y), min(a.z, b.z) };
    //});
    //if (Float3Length(temp_min - m_cell_box_min) > 0.001f)
    //{
    //    throw std::exception("Incorrect box!");
    //}
    //points_view.synchronize();
    //float_3 temp_max = amp_algorithms::reduce(points_view,
    //    [](const float_3& a, const float_3& b) restrict(cpu, amp) -> float_3
    //{
    //    return{ max(a.x, b.x), max(a.y, b.y), max(a.z, b.z) };
    //});
    //if (Float3Length(temp_max - m_cell_box_max) > 0.001f)
    //{
    //    throw std::exception("Incorrect box!");
    //}
    //m_cell_box_min = temp_min;
    //m_cell_box_max = temp_max;
}

void ParticleSystem::CreateCells()
{
    std::vector<float_3> a = *m_points;
    m_cell_size = int_3((m_cell_box_max - m_cell_box_min) / KERNEL);
    m_cell_size++;
    int cellNumber = m_cell_size.x * m_cell_size.y * m_cell_size.z;
    m_cells = std::make_unique<array<int>>(cellNumber, m_accelerator_view);
    array<int>& cells = *m_cells;
    array<float_3>& points = *m_points;
    array<int>& points_next = *m_points_next;
    parallel_for_each(m_cells->extent, [&cells](index<1> idx) restrict(amp) 
    {
        cells[idx] = -1;
    });
    float_3 box_min = this->m_cell_box_min;
    float_3 box_max = this->m_cell_box_max;
    int_3 cell_size = this->m_cell_size;
    float kernel = this->KERNEL;
    parallel_for_each(m_points->extent, [=, &points, &points_next, &cells](index<1> idx) restrict(amp)
    {
        int_3 cellIndex((points[idx] - box_min) / kernel);
        int i = (cellIndex.x * cell_size.y + cellIndex.y) * cell_size.z + cellIndex.z;
        points_next[idx] = Concurrency::atomic_exchange(&cells[i], idx[0]);
    });
}

void ParticleSystem::TestNeighbors()
{
    std::vector<int> indice = *m_neighbor_indice;
    std::vector<int> neighbor = *m_neighbors;
    std::vector<float_3> points = *m_points;
    std::vector<int> actual_indice;
    std::vector<int> actual_neighbors;
    actual_indice.push_back(0);
    for (int i = 0; i < points.size(); ++i)
    {
        for (int j = 0; j < points.size(); ++j)
        {
            if (Float3Length(points[j] - points[i]) < KERNEL)
            {
                actual_neighbors.push_back(j);
            }
        }
        actual_indice.push_back(actual_neighbors.size());
    }
    for (int i = 0; i < actual_indice.size(); ++i)
    {
        if (actual_indice[i] != indice[i])
            throw std::exception("Incorrect neighbors!");
    }
}

void ParticleSystem::TestDensity()
{
    std::vector<int> indice = *m_neighbor_indice;
    std::vector<int> neighbor = *m_neighbors;
    std::vector<float_3> points = *m_points;
    std::vector<float> density_list = *m_density;
    std::vector<float> actural_density_list;
    for (int i = 0; i < points.size(); ++i)
    {
        float density = 0;
        float_3 pi = points[i];
        for (int k = indice[i]; k < indice[i + 1]; ++k)
        {
            int j = neighbor[k];
            if (i == j) continue;
            float_3 pj = points[j];
            float_3 pij = pi - pj;
            float r2 = Float3LengthSq(pij);
            density += m_mass * m_poly6Value * pow(m_kernel2 - r2, 3);
        }
        density += m_selfDensity;
        actural_density_list.push_back(density);
    }
    for (int i = 0; i < density_list.size(); ++i)
    {
        if (abs(actural_density_list[i] - density_list[i]) > 0.01f)
            throw std::exception("Incorrect density!");
    }
}

void ParticleSystem::TestForce()
{
    static int time = 0;
    time++;
    std::vector<int> indice = *m_neighbor_indice;
    std::vector<int> neighbor = *m_neighbors;
    std::vector<float_3> points = *m_points;
    std::vector<float_3> velocity = *m_velocity_eval;
    std::vector<float> density_list = *m_density;
    std::vector<float_3> force = *m_force;
    std::vector<float> pressure = *m_pressure;
    std::vector<float_3> actural_force(force.size(), { 0, 0, 0 });
    std::vector<float> actural_pressure_list;
    for (auto p : density_list)
    {
        actural_pressure_list.push_back((pow(p / m_restDensity, 7) - 1) * m_gasConstant);
    }
    
    for (int i = 0; i < points.size(); ++i)
    {
        float_3 acceleration = { 0, 0, 0 };
        float_3 pi = points[i];
        float_3 vi = velocity[i];
        for (int k = indice[i]; k < indice[i + 1]; ++k)
        {
            int j = neighbor[k];
            if (i == j) continue;
            float_3 pj = points[j];
            float_3 pij = pi - pj;
            float r2 = Float3LengthSq(pij);
            float r = sqrt(r2);
            if (r < 0.001f)
                r = 0.001f;
            float V = m_mass / density_list[j] / 2;
            float kernel_r = m_kernel - r;
            float pressure_kernel = m_spikyValue * kernel_r * kernel_r;
            float temp_force = V * (pressure[i] + pressure[j]) * pressure_kernel;
            acceleration -= pij * temp_force / r;

            //float_3 vj = velocity[j];
            //float_3 vij = vi - vj;
            //float visc_kernel = m_viscoValue * kernel_r;
            //temp_force = V * m_viscosity * visc_kernel;
            //acceleration -= vij * temp_force;
        }
        actural_force[i] = acceleration;
    }
    for (int i = 0; i < points.size(); ++i)
    {
        if (Float3Length(actural_force[i] - force[i]) > 0.1f)
        {
            throw std::exception("Incorrect force!");
        }
    }
}

void ParticleSystem::GenerateNeighborTable()
{
    GetBox();
    CreateCells();
    array<int>& neighbor_indice(*m_neighbor_indice);
    array<float_3>& points(*m_points);
    array<int>& cells(*m_cells);
    array<int>& points_next(*m_points_next);
    neighbor_indice[m_points_vector.size()] = 1;
    int_3 cell_size = this->m_cell_size;
    float_3 box_min = this->m_cell_box_min;
    float kernel = this->KERNEL;
    parallel_for_each(points.extent,
        [=, &cells, &points, &neighbor_indice, &points_next](index<1> idx) restrict(amp) {
        int_3 cellIndex((points[idx] - box_min) / kernel);
        using Fluid::max;
        using Fluid::min;
        neighbor_indice[idx] = 0;
        if (idx[0] == 0) neighbor_indice[0] = 0;
        for (int x = max(cellIndex.x - 1, 0); x <= min(cellIndex.x + 1, cell_size.x - 1); ++x)
        {
            for (int y = max(cellIndex.y - 1, 0); y <= min(cellIndex.y + 1, cell_size.y - 1); ++y)
            {
                for (int z = max(cellIndex.z - 1, 0); z <= min(cellIndex.z + 1, cell_size.z - 1); ++z)
                {
                    int pos = (x * cell_size.y + y) * cell_size.z + z;
                    int i = cells[pos];
                    while (i != -1)
                    {
                        if (Float3Length(points[idx] - points[i]) < kernel)
                        {
                            Concurrency::atomic_fetch_inc(&neighbor_indice[idx]);
                        }
                        i = points_next[i];
                    }
                }
            }
        }
    });

    //std::vector<int> b = neighbor_indice;
    //Fluid::scan(neighbor_indice);
    //auto view = Concurrency::array_view<int>(neighbor_indice);
    //auto result = array<int>(neighbor_indice.extent);
    amp_algorithms::scan_exclusive(m_accelerator_view,
        neighbor_indice.view_as(neighbor_indice.extent),
        neighbor_indice.view_as(neighbor_indice.extent));
    //view.synchronize();
    //std::vector<int> a = result;
    int size = neighbor_indice[m_points_vector.size()];
    m_neighbors = std::make_unique<array<int>>(size, m_accelerator_view);
    array<int>& neighbors(*m_neighbors);
    parallel_for_each(points.extent,
        [=, &cells, &points, &neighbor_indice, &points_next, &neighbors](index<1> idx) restrict(amp) {
        int_3 cellIndex((points[idx] - box_min) / kernel);
        using Fluid::max;
        using Fluid::min;
        int j = neighbor_indice[idx];
        for (int x = max(cellIndex.x - 1, 0); x <= min(cellIndex.x + 1, cell_size.x - 1); ++x)
        {
            for (int y = max(cellIndex.y - 1, 0); y <= min(cellIndex.y + 1, cell_size.y - 1); ++y)
            {
                for (int z = max(cellIndex.z - 1, 0); z <= min(cellIndex.z + 1, cell_size.z - 1); ++z)
                {
                    int pos = (x * cell_size.y + y) * cell_size.z + z;
                    int i = cells[pos];
                    while (i != -1)
                    {
                        if (Float3Length(points[idx] - points[i]) < kernel)
                        {
                            neighbors[j++] = i;
                        }
                        i = points_next[i];
                    }
                }
            }
        }
    });

}

void ParticleSystem::ComputeDensity()
{
    parallel_for_each(m_points->extent,
        DensityFunction<Poly6Kernel>(KERNEL, this)
        );
}

template<typename KernelFunction>
ParticleSystem::DensityFunction<KernelFunction>::DensityFunction(float kernel, ParticleSystem* system):
    KernelFunction(kernel),
    m_points(system->Points()),
    m_neighbors(system->Neighbors()),
    m_neighborIndices(system->NeighborIndice()),
    m_density(system->Density())
{
}

template<typename KernelFunction>
void ParticleSystem::DensityFunction<KernelFunction>::operator() (index<1> idx) const restrict(amp)
{
    int i = idx[0];
    m_density[idx] = 0;
    for (int k = m_neighborIndices[i]; k < m_neighborIndices[i + 1]; ++k)
    {
        int j = m_neighbors[k];
        m_density[idx] += m_mass * Calculate(m_points[i], m_points[j]);
    }
}

template<typename KernelFunction>
ParticleSystem::PressureFunction<KernelFunction>::PressureFunction(float kernel, ParticleSystem* system) :
    KernelFunction(kernel),
    m_points(system->Points()), m_neighbors(system->Neighbors()),
    m_neighborIndices(system->NeighborIndice()), m_density(system->Density()),
    m_force(system->Force()), m_pressure(system->Pressure())
{
    
}

template<typename KernelFunction>
void ParticleSystem::PressureFunction<KernelFunction>::operator() (index<1> idx) const restrict(amp)
{
    int i = idx[0];
    float_3 ri = m_points[i];
    m_force[idx] = { 0, 0, 0 };
    for (int k = m_neighborIndices[i]; k < m_neighborIndices[i + 1]; ++k)
    {
        int j = m_neighbors[k];
        if (i == j) continue;
        float_3 rj = m_points[j];
        float V = m_mass / m_density[j] / 2;
        m_force[idx] -= V * (m_pressure[i] + m_pressure[j]) * CalculateGrad(ri, rj);
    }
}

template<typename KernelFunction>
ParticleSystem::ViscosityFunction<KernelFunction>::ViscosityFunction(float kernel, ParticleSystem* system) :
    KernelFunction(kernel),
    m_points(system->Points()), m_neighbors(system->Neighbors()),
    m_neighborIndices(system->NeighborIndice()), m_density(system->Density()),
    m_force(system->Force()), m_velocity_eval(system->VelocityEval())
{

}

template<typename KernelFunction>
void ParticleSystem::ViscosityFunction<KernelFunction>::operator() (index<1> idx) const restrict(amp)
{
    int i = idx[0];
    float_3 ri = m_points[i];
    float_3 vi = m_velocity_eval[i];
    for (int k = m_neighborIndices[i]; k < m_neighborIndices[i + 1]; ++k)
    {
        int j = m_neighbors[k];
        if (i == j) continue;
        float_3 rj = m_points[j];
        float_3 vj = m_velocity_eval[j];
        float V = m_mass / m_density[j] / 2;
        m_force[i] += (vj - vi) * V * 9.5f * CalculateLplc(ri, rj);
    }
}

//template<typename KernelFunction>
//ParticleSystem::ColorGradFunction<KernelFunction>::ColorGradFunction(ParticleSystem* system) :
//    KernelFunction(system->Kernel()), m_neighbors(system->Neighbors()),
//    m_neighbor_indices(system->NeighborIndice()), m_density(system->Density()),
//    m_color_grad(system->ColorGrad()), m_points(system->Points())
//{
//
//}
//
//template<typename KernelFunction>
//void ParticleSystem::ColorGradFunction<KernelFunction>::operator() (index<1> idx) const restrict(amp)
//{
//    int i = idx[0];
//    float_3 ri = m_points[i];
//    m_color_grad[i] = { 0, 0, 0 };
//    for (int k = m_neighbor_indices[i]; k < m_neighbor_indices[i + 1]; ++k)
//    {
//        int j = m_neighbors[k];
//        if (i == j) continue;
//        float_3 rj = m_points[j];
//        m_color_grad[idx] += MASS / density[j] * CalculateGrad(ri, rj);
//    }
//}

template<typename KernelFunction>
ParticleSystem::SurfaceTensionFunction<KernelFunction>::SurfaceTensionFunction(ParticleSystem* system):
    KernelFunction(system->Kernel()), m_neighbors(system->Neighbors()),
    m_neighbor_indices(system->NeighborIndice()), m_density(system->Density()),
    m_force(system->Force()), m_points(system->Points()), MASS(system->ParticleMass())
{

}

template<typename KernelFunction>
void ParticleSystem::SurfaceTensionFunction<KernelFunction>::operator() (index<1> idx) const restrict(amp)
{
    int i = idx[0];
    float_3 ri = m_points[i];
    float_3 color_grad = { 0, 0, 0 };
    float color_lplc = 0;

    for (int k = m_neighbor_indices[i]; k < m_neighbor_indices[i + 1]; ++k)
    {
        int j = m_neighbors[k];
        if (i == j) continue;
        float_3 rj = m_points[j];
        color_grad += MASS / m_density[j] * CalculateGrad(ri, rj);
        color_lplc += MASS / m_density[j] * CalculateLplc(ri, rj);
    }
    m_force[idx] += -SIGMA * color_lplc * color_grad / Float3Length(color_grad);
}

void ParticleSystem::ComputePressure()
{
    array<float>& pressure = *m_pressure;
    array<float>& density = *m_density;
    float gas_constant = m_gasConstant;
    float rest_density = m_restDensity;
    //std::vector<float> b = density;
    parallel_for_each(m_points->extent,
        [&pressure, &density, gas_constant, rest_density](index<1> idx) restrict(amp)
    {
        pressure[idx] = Concurrency::fast_math::powf(density[idx] / rest_density, 7) - 1;
        pressure[idx] *= gas_constant;
    });
    //std::vector<float> a = pressure;
    //    p.pressure = (pow(p.density / m_restDensity, 7) - 1) * m_gasConstant;
    parallel_for_each(m_points->extent,
        PressureFunction<SpikyKernel>(KERNEL, this));
}

void ParticleSystem::ComputeViscosity()
{
    parallel_for_each(m_points->extent,
        ViscosityFunction<ViscosityKernel>(KERNEL, this));
}

void ParticleSystem::ComputeSurfaceTension()
{
    parallel_for_each(m_points->extent,
        SurfaceTensionFunction<Poly6Kernel>(this));
}

void ParticleSystem::ComputeInternalForce()
{
    //using DirectX::XMVECTOR;
    //for (Particle& p : m_particles)
    //{
    //    p.pressure = (pow(p.density / m_restDensity, 7) - 1) * m_gasConstant;
    //}

    //for (Particle& p : m_particles)
    //{
    //    XMVECTOR pi = DirectX::XMLoadFloat3(&p.position);
    //    XMVECTOR vi = DirectX::XMLoadFloat3(&p.velocityEval);
    //    XMVECTOR acceleration = DirectX::XMLoadFloat3(&p.acceleration);
    //    for (int j : p.neighbors)
    //    {
    //        XMVECTOR pj = DirectX::XMLoadFloat3(&m_particles[j].position);
    //        XMVECTOR pij = DirectX::XMVectorSubtract(pi, pj);
    //        float r2;
    //        DirectX::XMStoreFloat(&r2, DirectX::XMVector3LengthSq(pij));
    //        float r = sqrt(r2);
    //        if (r < 0.001f) r = 0.001f;
    //        float V = m_mass / m_particles[j].density / 2;
    //        float kernel_r = m_kernel - r;
    //        float pressureKernel = m_spikyValue * kernel_r * kernel_r;
    //        float tempForce = V * (p.pressure + m_particles[j].pressure) * pressureKernel;
    //        acceleration = DirectX::XMVectorSubtract(acceleration,
    //            DirectX::XMVectorScale(pij, tempForce / r));

    //        XMVECTOR vj = DirectX::XMLoadFloat3(&m_particles[j].velocityEval);
    //        XMVECTOR vij = DirectX::XMVectorSubtract(vi, vj);
    //        float viscKernel = m_viscoValue * kernel_r;
    //        tempForce = V * m_viscosity * viscKernel;
    //        acceleration = DirectX::XMVectorSubtract(acceleration,
    //            DirectX::XMVectorScale(vij, tempForce));
    //    }
    //    DirectX::XMStoreFloat3(&p.acceleration, acceleration);
    //}
    array<float_3>& force = *m_force;
    parallel_for_each(force.extent,
        [&force](index<1> idx) restrict(amp)
    {
        force[idx] = { 0, 0, 0 };
    });
    ComputePressure();
    ComputeViscosity();
    ComputeSurfaceTension();
    //std::vector<float_3> f = *m_force;
}

void ParticleSystem::Advect()
{
    float BOUNDARY = 0.0001f;
    array<float_3>& points = *m_points;
    array<float>& density = *m_density;
    array<float_3>& force = *m_force;
    array<float_3>& velocity = *m_velocity;
    array<float_3>& velocity_eval = *m_velocity_eval;
    float time_step = m_timeStep;
    float wall_damping = m_wallDamping;
    float_3 box_max = m_box_max;
    //std::vector<float_3> f = force;
    //std::vector<float> den = density;
    parallel_for_each(m_points->extent,
        [=, &points, &density, &force, &velocity, &velocity_eval] (index<1> idx) restrict(amp)
    {
        float_3 acceleration = force[idx] / density[idx];
        if (Float3Length(acceleration) > 3000.0f)
        {
            acceleration *= 3000.0f / Float3Length(acceleration);
        }
        acceleration += { 0, -9.8f, 0 };
        velocity_eval[idx] = velocity[idx];
        velocity[idx] += acceleration * time_step;
        points[idx] += velocity[idx] * time_step;
        if (points[idx].x >= box_max.x)
        {
            velocity[idx].x *= wall_damping;
            points[idx].x = box_max.x;
        }
        if (points[idx].y >= box_max.y)
        {
            velocity[idx].y *= wall_damping;
            points[idx].y = box_max.y;
        }
        if (points[idx].z >= box_max.z)
        {
            velocity[idx].z *= wall_damping;
            points[idx].z = box_max.z;
        }
        if (points[idx].x < 0.0f)
        {
            velocity[idx].x *= wall_damping;
            points[idx].x = 0.0f;
        }
        if (points[idx].y < 0.0f)
        {
            velocity[idx].y *= wall_damping;
            points[idx].y = 0.0f;
        }
        if (points[idx].z < 0.0f)
        {
            velocity[idx].z *= wall_damping;
            points[idx].z = 0.0f;
        }
        velocity_eval[idx] += velocity[idx];
        velocity_eval[idx] /= 2;
    });
    //std::vector<float_3> vel = velocity;
}
