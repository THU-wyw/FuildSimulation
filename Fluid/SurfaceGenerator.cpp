#include "pch.h"
#include "SurfaceGenerator.h"
#include "Content\MarchingCubesRes.h"
#include <fstream>
using namespace Fluid;
using Fluid::ParticleSystem;
using Concurrency::graphics::float_3;
using Concurrency::array;
using Concurrency::index;
using Concurrency::graphics::int_3;
using namespace DirectX;

SurfaceGenerator::SurfaceGenerator()
{
}

void SurfaceGenerator::Generate(std::vector<Vertex::Basic32> &vertices)
{
	vertices.clear();
	// m_system->GenerateNeighborTable();
	XMFLOAT3 cell_box_min = ConvertToDX(m_system->m_cell_box_min);
	// std::fstream fout("output.txt", std::ios::out);
	std::vector<Concurrency::graphics::float_3> list = m_system->Points();
	std::vector<XMFLOAT3> points(list.size());
	std::transform(list.begin(),
		list.end(),
		points.begin(),
		[](const Concurrency::graphics::float_3& p) -> XMFLOAT3 {
		return{ p.x, p.y, p.z };
	});
	int pn = points.size();
	int grid_n = 2;
	float grid_step = m_system->Kernel() / (float)grid_n;
	XMINT3 cell_size = ConvertToDX(m_system->m_cell_size);
	// padding the cells
	cell_size.x = (cell_size.x + 2) * grid_n;
	cell_size.y = (cell_size.y + 2) * grid_n;
	cell_size.z = (cell_size.z + 2) * grid_n;
	XMINT3 grid_size(cell_size.x + 1, cell_size.y + 1, cell_size.z + 1);
	int cell_num = cell_size.x * cell_size.y * cell_size.z;
	int grid_num = grid_size.x * grid_size.y * grid_size.z;
	std::vector<float> grid_density(grid_num);
	std::vector<XMFLOAT3> grid_normal(grid_num);
	// calculate cells and points_next
	std::vector<int> cells(cell_num, -1);
	std::vector<int> points_next(points.size(), -1);
	for (int i = 0; i < pn; i++) {
		int x = int((points[i].x - cell_box_min.x) / grid_step) + grid_n;
		if (x < 0 || x >= cell_size.x) continue;
		// x = x >= cell_size.x ? cell_size.x - 1 : x;
		int y = int((points[i].y - cell_box_min.y) / grid_step) + grid_n;
		if (y < 0 || y >= cell_size.y) continue;
		// y = y >= cell_size.y ? cell_size.y - 1 : y;
		int z = int((points[i].z - cell_box_min.z) / grid_step) + grid_n;
		if (z < 0 || z >= cell_size.z) continue;
		// z = z >= cell_size.z ? cell_size.z - 1 : z;
		int cellIdx = GridIndex(cell_size, x, y, z);
		if (cells[cellIdx] < 0) {
			cells[cellIdx] = i;
		}
		else {
			int j;
			for (j = cells[cellIdx]; points_next[j] >= 0; j = points_next[j]);
			points_next[j] = i;
		}
	}
	// calculate grid density
	for (int i = 0; i < grid_num; i++) {
		int index = i;
		int z = index % grid_size.z; index /= grid_size.z;
		int y = index % grid_size.y; index /= grid_size.y;
		int x = index % grid_size.x;
		float density = 0.0f;
		XMFLOAT3 grid_pos((x - 1) * grid_step + cell_box_min.x, (y - 1) * grid_step + cell_box_min.y, (z - 1) * grid_step + cell_box_min.z);
		for (int ix = std::max(0, x - grid_n); ix <= std::min(x + grid_n - 1, cell_size.x - 1); ix++) {
			for (int iy = std::max(0, y - grid_n); iy <= std::min(y + grid_n - 1, cell_size.y - 1); iy++) {
				for (int iz = std::max(0, z - grid_n); iz <= std::min(z + grid_n - 1, cell_size.z - 1); iz++) {
					int cellIdx = GridIndex(cell_size, ix, iy, iz);
					for (int j = cells[cellIdx]; j >= 0; j = points_next[j]) {
						double d = kernel->Calculate(grid_pos, points[j]);
						if (d > 0) density += d;
					}
				}
			}
		}
		grid_density[i] = density;
	}
	// calculate grid normal
	for (int i = 0; i < grid_num; i++) {
		int index = i;
		int z = index % grid_size.z; index /= grid_size.z;
		int y = index % grid_size.y; index /= grid_size.y;
		int x = index % grid_size.x;
		float xl = x == 0 ? 0 : grid_density[GridIndex(grid_size, x - 1, y, z)];
		float xr = x == grid_size.x - 1 ? 0 : grid_density[GridIndex(grid_size, x + 1, y, z)];
		float yl = y == 0 ? 0 : grid_density[GridIndex(grid_size, x, y - 1, z)];
		float yr = y == grid_size.y - 1 ? 0 : grid_density[GridIndex(grid_size, x, y + 1, z)];
		float zl = z == 0 ? 0 : grid_density[GridIndex(grid_size, x, y, z - 1)];
		float zr = z == grid_size.z - 1 ? 0 : grid_density[GridIndex(grid_size, x, y, z + 1)];
		XMVECTOR tmp = XMVectorSet(xl - xr, yl - yr, zl - zr, 0.0f);
		// tmp = XMVector3Normalize(tmp);
		XMStoreFloat3(&(grid_normal[i]), tmp);
	}
	float threshold = 10.0f;
	// calculate vertices
	for (int i = 0; i < cell_num; i++) {
		int index = i;
		int z = index % cell_size.z; index /= cell_size.z;
		int y = index % cell_size.y; index /= cell_size.y;
		int x = index % cell_size.x;
		// calculate vert value
		unsigned int verts = 0;
		for (int j = 0; j < 8; j++) {
			XMUINT3 dpos = cellPos[j];
			int d = GridIndex(grid_size, x + dpos.x, y + dpos.y, z + dpos.z);
			if (grid_density[d] < threshold) verts |= (1 << j);
		}
		XMVECTOR base_pos = XMVectorSet(x - 1, y - 1, z - 1, 0.0);
		for (int j = 0, n = numVertsTable[verts]; j < n; j++) {
			unsigned int edge_index = triTable[verts][j];
			XMUINT2 pair_index = idxTable[edge_index];
			Vertex::Basic32 vert;
			// point position
			XMUINT3 e1 = cellPos[pair_index.x], e2 = cellPos[pair_index.y];
			XMVECTOR ev1 = XMLoadUInt3(&e1), ev2 = XMLoadUInt3(&e2);
			XMVECTOR pos = (base_pos + (ev1 + ev2) / 2.0f) * grid_step;
			// normal
			int ni1 = GridIndex(grid_size, x + e1.x, y + e1.y, z + e1.z);
			int ni2 = GridIndex(grid_size, x + e2.x, y + e2.y, z + e2.z);
			XMVECTOR n1 = XMLoadFloat3(&grid_normal[ni1]);
			XMVECTOR n2 = XMLoadFloat3(&grid_normal[ni2]);
			XMVECTOR normal = (n1 + n2) / 2.0f;
			normal = XMVector3Normalize(normal);
			XMStoreFloat3(&vert.Pos, pos);
			XMStoreFloat3(&vert.Normal, normal);
			vertices.push_back(vert);
		}
	}
}