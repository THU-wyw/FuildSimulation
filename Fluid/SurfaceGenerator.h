#pragma once
#include "Sph\ParticleSystem.h"
#include "Sph\KernelFuntions.h"
#include "Content\ShaderStructures.h"
#include "External\CubeMap\Vertex.h"
namespace Fluid
{
	class SurfaceGenerator
	{
	public:
		SurfaceGenerator();
		void Init(ParticleSystem *system) { 
			m_system = system; 
			kernel = new Poly6Kernel(m_system->Kernel()); 
		}
		void Generate(std::vector<Vertex::Basic32> &vertices);
	private:
		inline DirectX::XMFLOAT3 ConvertToDX(Concurrency::graphics::float_3 f3){
			return DirectX::XMFLOAT3(f3.x, f3.y, f3.z);
		}
		inline DirectX::XMINT3 ConvertToDX(Concurrency::graphics::int_3 i3) {
			return DirectX::XMINT3(i3.x, i3.y, i3.z);
		}
		inline int GridIndex(DirectX::XMINT3 cell_size, int x, int y, int z) {
			return (x * cell_size.y + y) * cell_size.z + z;
		}
		Poly6Kernel *kernel;
		ParticleSystem *m_system;
		int x;
	};
};


