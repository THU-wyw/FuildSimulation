#include "pch.h"
#include "Sample3DSceneRenderer.h"
#include "..\Common\DirectXHelper.h"
#include <vector>
#include <DDSTextureLoader.h>
#include <fstream>

using namespace Fluid;
using namespace DirectX;
using namespace Windows::Foundation;

void Sample3DSceneRenderer::InitEnvMapShaders()
{
	// Load shaders asynchronously.
	auto loadVSTask = DX::ReadDataAsync(L"EnvMapVertexShader.cso");
	auto loadPSTask = DX::ReadDataAsync(L"EnvMapPixelShader.cso");

	// After the vertex shader file is loaded, create the shader and input layout.
	auto createVSTask = loadVSTask.then([this](const std::vector<byte>& fileData) {
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateVertexShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_envmap_vertexShader
				)
			);
		static const D3D11_INPUT_ELEMENT_DESC vertexDesc[] =
		{
			{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		};

		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateInputLayout(
				vertexDesc,
				ARRAYSIZE(vertexDesc),
				&fileData[0],
				fileData.size(),
				&m_envmap_layout
				)
			);
	});

	// After the pixel shader file is loaded, create the shader and constant buffer.
	auto createPSTask = loadPSTask.then([this](const std::vector<byte>& fileData) {
		ID3D11Device* mDevice = m_deviceResources->GetD3DDevice();
		
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreatePixelShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_envmap_pixelShader
				)
			);

		// setup constant buffer
		D3D11_BUFFER_DESC Desc;
		Desc.Usage = D3D11_USAGE_DYNAMIC;
		Desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		Desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
		Desc.MiscFlags = 0;
		Desc.ByteWidth = sizeof(cbPerFrame2);
		DX::ThrowIfFailed( 
			mDevice->CreateBuffer(&Desc, NULL, &m_envmap_vertexCB)
		);

		// setup sampler state 
		D3D11_SAMPLER_DESC SamDesc;
		SamDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
		SamDesc.AddressU = D3D11_TEXTURE_ADDRESS_WRAP;
		SamDesc.AddressV = D3D11_TEXTURE_ADDRESS_WRAP;
		SamDesc.AddressW = D3D11_TEXTURE_ADDRESS_WRAP;
		SamDesc.MipLODBias = 0.0f;
		SamDesc.MaxAnisotropy = 1;
		SamDesc.ComparisonFunc = D3D11_COMPARISON_ALWAYS;
		SamDesc.BorderColor[0] = SamDesc.BorderColor[1] = SamDesc.BorderColor[2] = SamDesc.BorderColor[3] = 0;
		SamDesc.MinLOD = 0;
		SamDesc.MaxLOD = D3D11_FLOAT32_MAX;
		DX::ThrowIfFailed(
			mDevice->CreateSamplerState(&SamDesc, &m_envmap_samplerState)
		);

		// get environment map cube texture
		DX::ThrowIfFailed(
			DirectX::CreateDDSTextureFromFile(mDevice, L"Media\\uffizi_cross.dds", NULL, &m_envMapSRV)
		);

		// depth stencil state
		D3D11_DEPTH_STENCIL_DESC DSDesc;
		ZeroMemory(&DSDesc, sizeof(D3D11_DEPTH_STENCIL_DESC));
		DSDesc.DepthEnable = TRUE;
		DSDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
		DSDesc.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
		DSDesc.StencilEnable = FALSE;
		DX::ThrowIfFailed(
			mDevice->CreateDepthStencilState(&DSDesc, &m_envmap_depthSS)
		);
	});

	(createPSTask && createVSTask).then([this]() {
		m_envmapComplete = true;
	});
}

void Sample3DSceneRenderer::RenderEnvMap()
{
	if (!m_modelComplete)
	{
		return;
	}
	auto context = m_deviceResources->GetD3DDeviceContext();
	ID3D11Device* mDevice = m_deviceResources->GetD3DDevice();

	// map texels to pixels
	Size wSize = m_deviceResources->GetLogicalSize();
	float fHighW = -1.0f - (1.0f / (float)wSize.Width);
	float fHighH = -1.0f - (1.0f / (float)wSize.Height);
	float fLowW = 1.0f + (1.0f / (float)wSize.Width);
	float fLowH = 1.0f + (1.0f / (float)wSize.Height);
	XMFLOAT4* pVertices = new XMFLOAT4[6];

	pVertices[0] = XMFLOAT4(fLowW, fLowH, 0.0f, 1.0f);
	pVertices[1] = XMFLOAT4(fLowW, fHighH, 0.0f, 1.0f);
	pVertices[2] = XMFLOAT4(fHighW, fLowH, 0.0f, 1.0f);
	pVertices[3] = XMFLOAT4(fHighW, fHighH, 0.0f, 1.0f);
	
	UINT uiVertBufSize = 4 * sizeof(XMFLOAT4);
	//Vertex Buffer
	D3D11_BUFFER_DESC vbdesc;
	vbdesc.ByteWidth = uiVertBufSize;
	vbdesc.Usage = D3D11_USAGE_IMMUTABLE;
	vbdesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vbdesc.CPUAccessFlags = 0;
	vbdesc.MiscFlags = 0;

	D3D11_SUBRESOURCE_DATA InitData;
	InitData.pSysMem = pVertices;
	mDevice->CreateBuffer(&vbdesc, &InitData, &m_pVB);
	if (pVertices)
	{
		delete[] pVertices;
	}
	pVertices = nullptr;

	UINT uStrides = sizeof(XMFLOAT4);
	UINT uOffsets = 0;
	context->IASetVertexBuffers(0, 1, m_pVB.GetAddressOf(), &uStrides, &uOffsets);
	context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
	//context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	context->IASetInputLayout(m_envmap_layout.Get());
	context->VSSetShader(m_envmap_vertexShader.Get(), NULL, 0);
	context->PSSetShader(m_envmap_pixelShader.Get(), NULL, 0);

	D3D11_MAPPED_SUBRESOURCE MappedResource;
	context->Map(m_envmap_vertexCB.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource);
	cbPerFrame2* pVSPerObject = (cbPerFrame2*)MappedResource.pData;
	//pVSPerObject->mWordViewProj = XMMatrixTranspose(XMMatrixInverse(NULL, m_vertexConstants.WVP));
	pVSPerObject->mWordViewProj = m_vertexConstants.WVP;

	context->Unmap(m_envmap_vertexCB.Get(), 0);
	context->VSSetConstantBuffers(0, 1, m_envmap_vertexCB.GetAddressOf());
	context->PSSetSamplers(0, 1, m_envmap_samplerState.GetAddressOf());
	context->PSSetShaderResources(0, 1, m_envMapSRV.GetAddressOf());
	context->OMSetDepthStencilState(m_envmap_depthSS.Get(), 0);
	context->Draw(4, 0);
}