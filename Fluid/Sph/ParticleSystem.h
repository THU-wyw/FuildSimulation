#pragma once
#ifndef SPH_PARTICLE_SYSTEM_H
#define SPH_PARTICLE_SYSTEM_H

#include <vector>
#include <DirectXCollision.h>
#include <d3d11_1.h>
#include <amp_graphics.h>
#include <amp.h>

#include "ParticleGrid.h"

namespace Fluid
{
	class SurfaceGenerator;
    class ParticleSystem
    {
    public:
        ParticleSystem(ID3D11Device2* device);
        ~ParticleSystem();
        void InitParticles();
        void Update();
        HRESULT GetDataBuffer(void** d3dbuffer) const
        {
            return Concurrency::direct3d::get_buffer(*m_points)->QueryInterface(__uuidof(ID3D11Buffer), (LPVOID*)d3dbuffer);
        }
        const std::vector<Particle>& GetParticles() const { return m_particles; }
        inline float Kernel() const{ return KERNEL; }
        inline float ParticleMass() const { return m_mass; }
        inline Concurrency::array<Concurrency::graphics::float_3>& Points() { return *m_points; }
        inline Concurrency::array<Concurrency::graphics::float_3>& Velocity() { return *m_velocity; }
        inline Concurrency::array<Concurrency::graphics::float_3>& VelocityEval() { return *m_velocity_eval; }
        inline Concurrency::array<Concurrency::graphics::float_3>& Force() { return *m_force; }
        inline Concurrency::array<float>& Density() { return *m_density; }
        inline Concurrency::array<float>& Pressure() { return *m_pressure; }
        inline Concurrency::array<int>& Neighbors() { return *m_neighbors; }
        inline Concurrency::array<int>& NeighborIndice() { return *m_neighbor_indice; }
        inline Concurrency::array<Concurrency::graphics::float_3>& ColorGrad() { return *m_color_grad; }
        inline Concurrency::array<float>& ColorLplc() { return *m_color_lplc; }

		friend class SurfaceGenerator;
    private:
        template<typename KernelFunction>
        class DensityFunction : public KernelFunction
        {
        public:
            DensityFunction(
                float kernel,
                ParticleSystem* system);
            DensityFunction(const DensityFunction&) = delete;
            DensityFunction(const DensityFunction&&) = delete;
            void operator() (Concurrency::index<1> idx) const restrict(amp);
        private:
            float m_mass = 0.02f;
            const Concurrency::array<Concurrency::graphics::float_3>& m_points;
            const Concurrency::array<int>& m_neighbors;
            const Concurrency::array<int>& m_neighborIndices;
            Concurrency::array<float>& m_density;
        };

        template<typename KernelFunction>
        class PressureFunction : public KernelFunction
        {
        public:
            PressureFunction(
                float kernel,
                ParticleSystem* system);
            void operator() (Concurrency::index<1> idx) const restrict(amp);
        private:
            float m_mass = 0.02f;
            const Concurrency::array<Concurrency::graphics::float_3>& m_points;
            const Concurrency::array<int>& m_neighbors;
            const Concurrency::array<int>& m_neighborIndices;
            const Concurrency::array<float>& m_pressure;
            const Concurrency::array<float>& m_density;
            Concurrency::array <Concurrency::graphics::float_3>& m_force;
        };

        template<typename KernelFunction>
        class ViscosityFunction : KernelFunction
        {
        public:
            ViscosityFunction(
                float kernel,
                ParticleSystem* system);
            void operator() (Concurrency::index<1> idx) const restrict(amp);
        private:
            float m_mass = 0.02f;
            const float VISCOSITY = 9.5f;
            const Concurrency::array<Concurrency::graphics::float_3>& m_points;
            const Concurrency::array<Concurrency::graphics::float_3>& m_velocity_eval;
            const Concurrency::array<int>& m_neighbors;
            const Concurrency::array<int>& m_neighborIndices;
            const Concurrency::array<float>& m_density;
            Concurrency::array<Concurrency::graphics::float_3>& m_force;
        };

        //template<typename KernelFunction>
        //class ColorGradFunction : public KernelFunction
        //{
        //public:
        //    ColorGradFunction(ParticleSystem* system);
        //    void operator() (Concurrency::index<1> idx) const restrict(amp);
        //private:
        //    const float MASS;
        //    const Concurrency::array<Concurrency::graphics::float_3>& m_points;
        //    const Concurrency::array<int>& m_neighbors;
        //    const Concurrency::array<int>& m_neighbor_indices;
        //    const Concurrency::array<float>& m_density;
        //    Concurrency::array<float_3>& m_color_grad;
        //};

        //template<typename KernelFunction>
        //class ColorLplcFunction : public KernelFunction
        //{
        //public:
        //    ColorLplcFunction(ParticleSystem* system);
        //    void operator() (Concurrency::index<1> idx) const restrict(amp);
        //private:
        //    const float MASS;
        //    const Concurrency::array<Concurrency::graphics::float_3>& m_points;
        //    const Concurrency::array<int>& m_neighbors;
        //    const Concurrency::array<int>& m_neighbor_indices;
        //    const Concurrency::array<float>& m_density;
        //    Concurrency::array<float>& m_color_lplc;
        //};

        template<typename KernelFunction>
        class SurfaceTensionFunction : public KernelFunction
        {
        public:
            SurfaceTensionFunction(ParticleSystem* system);
            void operator() (Concurrency::index<1> idx) const restrict(amp);

        private:
            float const SIGMA = 1.0f;
            const float MASS;
            const Concurrency::array<Concurrency::graphics::float_3>& m_points;
            const Concurrency::array<int>& m_neighbors;
            const Concurrency::array<int>& m_neighbor_indices;
            const Concurrency::array<float>& m_density;
            Concurrency::array<Concurrency::graphics::float_3>& m_force;
        };
        void GenerateNeighborTable();
        void ComputeInternalForce();
        void ComputeDensity();
        void ComputePressure();
        void ComputeViscosity();
        void ComputeSurfaceTension();
        void Advect();
        void GetBox();
        void CreateCells();
        void TestNeighbors();
        void TestDensity();
        void TestForce();
        std::vector<Particle> m_particles;
        DirectX::XMFLOAT3 m_boxMin;
        DirectX::XMFLOAT3 m_boxMax;
        const DirectX::XMFLOAT3 m_gravity = { 0, -9.8f, 0 };

        const float m_kernel = 0.04f;
        const float m_mass = 0.02f;
        const float m_wallDamping = -0.5f;
        const float m_restDensity = 1000.0f;
        const float m_gasConstant = 1.0f;
        const float m_viscosity = 9.5f;
        const float m_timeStep = 0.003f;
        const float m_poly6Value = 315.0f / (64.0f * DirectX::XM_PI * pow(m_kernel, 9));
        const float m_spikyValue = -45.0f / (DirectX::XM_PI * pow(m_kernel, 6));
        const float m_viscoValue = 45.0f / (DirectX::XM_PI * pow(m_kernel, 6));
        const float m_poly6Grad = -945.0f / (32 * DirectX::XM_PI * pow(m_kernel, 9));
        const float m_poly6Lplc = -945.0f / (8 * DirectX::XM_PI * pow(m_kernel, 9));
        const float m_kernel2 = m_kernel * m_kernel;
        const float m_selfDensity = m_mass * m_poly6Value * pow(m_kernel, 6);
        ParticleGrid m_grid;
        Concurrency::accelerator_view m_accelerator_view;
        const float KERNEL = 0.04f;
        const float MASS = 0.02f;
        Concurrency::graphics::int_3 m_cell_size;
        Concurrency::graphics::float_3 m_box_min;
        Concurrency::graphics::float_3 m_box_max;
        Concurrency::graphics::float_3 m_cell_box_min;
        Concurrency::graphics::float_3 m_cell_box_max;
        std::vector<Concurrency::graphics::float_3> m_points_vector;
        std::unique_ptr<Concurrency::array<Concurrency::graphics::float_3>> m_points;
        std::unique_ptr<Concurrency::array<Concurrency::graphics::float_3>> m_velocity;
        std::unique_ptr<Concurrency::array<Concurrency::graphics::float_3>> m_velocity_eval;
        std::unique_ptr<Concurrency::array<Concurrency::graphics::float_3>> m_force;
        std::unique_ptr<Concurrency::array<float>> m_density;
        std::unique_ptr<Concurrency::array<float>> m_pressure;
        std::unique_ptr<Concurrency::array<int>> m_cells;
        std::unique_ptr<Concurrency::array<int>> m_points_next;
        std::unique_ptr<Concurrency::array<int>> m_neighbors;
        std::unique_ptr<Concurrency::array<int>> m_neighbor_indice;
        std::unique_ptr<Concurrency::array<Concurrency::graphics::float_3>> m_color_grad;
        std::unique_ptr<Concurrency::array<float>> m_color_lplc;
    };
}

#endif