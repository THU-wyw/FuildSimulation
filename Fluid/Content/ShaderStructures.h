#pragma once

namespace Fluid
{
	// Constant buffer used to send MVP matrices to the vertex shader.
	struct ModelViewProjectionConstantBuffer
	{
		DirectX::XMFLOAT4X4 model;
		DirectX::XMFLOAT4X4 view;
		DirectX::XMFLOAT4X4 projection;
	};

	// Used to send per-vertex data to the vertex shader.
	struct VertexPositionColor
	{
		DirectX::XMFLOAT3 pos;
		DirectX::XMFLOAT3 color;
	};

	struct mVertex
	{
		DirectX::XMFLOAT3 pos;
		DirectX::XMFLOAT2 tex;
		DirectX::XMFLOAT3 normal;
	};

	struct VertexConstants
	{
		DirectX::XMMATRIX MVP;
		DirectX::XMMATRIX World;
	};

	struct Light
	{
		Light()
		{
			ZeroMemory(this, sizeof(Light));
		}
		DirectX::XMFLOAT3 pos;
		float range;
		DirectX::XMFLOAT3 dir;
		float cone;
		DirectX::XMFLOAT3 att;
		float pad2;
		DirectX::XMFLOAT4 ambient;
		DirectX::XMFLOAT4 diffuse;
	};

	struct cbPerFrame
	{
		Light  light;
	};
}