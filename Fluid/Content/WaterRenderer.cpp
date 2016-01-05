#include "pch.h"
#include "Sample3DSceneRenderer.h"
#include "..\Common\DirectXHelper.h"
#include <vector>
#include <DDSTextureLoader.h>
#include <fstream>

using namespace Fluid;

using namespace DirectX;
using namespace Windows::Foundation;

HRESULT CreateTypedTexture2D(ID3D11Device* pd3dDevice, 
	DXGI_FORMAT format, 
	DXGI_FORMAT formatUAV, 
	UINT iWidth, UINT iHeight, 
	ID3D11Texture2D** ppTex,
	ID3D11ShaderResourceView** ppSRV, 
	ID3D11RenderTargetView** ppRTV, 
	ID3D11UnorderedAccessView** ppUAV, 
	const void* pInitialData = NULL, 
	UINT iDepth = 1)
{
	HRESULT hr = S_OK;

	// Create SB
	if (ppTex) {
		D3D11_TEXTURE2D_DESC texDesc;
		ZeroMemory(&texDesc, sizeof(texDesc));
		texDesc.ArraySize = iDepth;
		texDesc.CPUAccessFlags = 0;
		texDesc.Format = format;
		texDesc.Height = iHeight;
		texDesc.MipLevels = 1;
		texDesc.SampleDesc.Count = 1;
		texDesc.SampleDesc.Quality = 0;
		texDesc.Width = iWidth;
		texDesc.Usage = D3D11_USAGE_DEFAULT;
		texDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET | (ppUAV ? D3D11_BIND_UNORDERED_ACCESS : 0);
		texDesc.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA texInitData;
		ZeroMemory(&texInitData, sizeof(texInitData));
		texInitData.pSysMem = pInitialData;
		DX::ThrowIfFailed(pd3dDevice->CreateTexture2D(&texDesc, (pInitialData) ? &texInitData : NULL, ppTex));
	}

	// Create SRV
	if (ppTex && ppSRV) {
		DX::ThrowIfFailed(pd3dDevice->CreateShaderResourceView(*ppTex, NULL, ppSRV));
	}

	// Create UAV
	if (ppUAV) {
		D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc;
		ZeroMemory(&uavDesc, sizeof(uavDesc));
		uavDesc.Format = formatUAV;
		uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
		uavDesc.Texture2D.MipSlice = 0;
		DX::ThrowIfFailed(pd3dDevice->CreateUnorderedAccessView(*ppTex, &uavDesc, ppUAV));
	}

	if (ppRTV) {
		DX::ThrowIfFailed(pd3dDevice->CreateRenderTargetView(*ppTex, NULL, ppRTV));
	}

	return hr;
}

void Sample3DSceneRenderer::RenderCaustic()
{
	ID3D11DeviceContext* mContext = m_deviceResources->GetD3DDeviceContext();
	mContext->OMSetRenderTargets(1, m_model_causticRTV.GetAddressOf(), NULL);
	FLOAT clrz[4] = { 0 };
	mContext->ClearRenderTargetView(m_model_causticRTV.Get(), clrz);

	// TO DO:
}

void Sample3DSceneRenderer::InitCausticTex()
{
	ID3D11Device* mDevice = m_deviceResources->GetD3DDevice();
	Size wSize = m_deviceResources->GetLogicalSize();
	float width = wSize.Width;
	float height = wSize.Height;

	DX::ThrowIfFailed(
		CreateTypedTexture2D(mDevice,
			DXGI_FORMAT_R16G16B16A16_FLOAT,		// format
			DXGI_FORMAT_R16G16B16A16_FLOAT,		// formatUAV
			width, // width
			height, // height
			m_model_caustic.GetAddressOf(), // ID3D11Texture2D** ppTex,
			m_model_causticSRV.GetAddressOf(), // ID3D11ShaderResourceView** ppSRV,
			m_model_causticRTV.GetAddressOf(), // ID3D11RenderTargetView** ppRTV,
			NULL // ID3D11UnorderedAccessView** ppUAV,
				 // const void* pInitialData = NULL,
		)
	);
}

// Called once per frame, rotates the cube and calculates the model and view matrices.
void Sample3DSceneRenderer::Update(DX::StepTimer const& timer)
{
	if (!m_loadingComplete || !m_modelComplete || !m_envmapComplete)
		return;
	if (!m_tracking)
	{
		// Convert degrees to radians, then convert seconds to rotation angle
		float radiansPerSecond = XMConvertToRadians(m_degreesPerSecond);
		double totalRotation = timer.GetTotalSeconds() * radiansPerSecond;
		float radians = static_cast<float>(fmod(totalRotation, XM_2PI));

		Rotate(0);
	}

	m_system.Update();

	std::vector<Concurrency::graphics::float_3> list = m_system.Points();
	std::transform(list.begin(),
		list.end(),
		m_points.begin(),
		[](const Concurrency::graphics::float_3& p) -> VertexPositionColor {
		return{ { p.x, p.y, p.z },{ 1.0f, 1.0f, 1.0f } };
	});
	m_surface.Generate(model_vertices);
	InitModelVertices();
	XMMATRIX world = XMMatrixScaling(1.0f, 1.0f, 1.0f);;
	world = XMMatrixMultiply(XMLoadFloat4x4(&m_constantBufferData.model), world);
	XMMATRIX view = (XMLoadFloat4x4(&m_constantBufferData.view));
	XMMATRIX proj = (XMLoadFloat4x4(&m_constantBufferData.projection));
	m_vertexConstants.WVP = XMMatrixMultiply(proj, XMMatrixMultiply(view, world));
	m_vertexConstants.World = world;
	m_vertexConstants.WV = XMMatrixMultiply(view, world);
}

void Sample3DSceneRenderer::InitModelShaders()
{
	constbuffPerFrame.light.pos = XMFLOAT3(0.0f, 0.7f, 1.5f);
	constbuffPerFrame.light.dir = XMFLOAT3(0.0f, -0.7f, -0.8f);
	constbuffPerFrame.light.range = 1000.0f;
	constbuffPerFrame.light.cone = 20.0f;
	constbuffPerFrame.light.att = XMFLOAT3(0.4f, 0.02f, 0.000f);
	constbuffPerFrame.light.ambient = XMFLOAT4(0.2f, 0.2f, 0.2f, 1.0f);
	constbuffPerFrame.light.diffuse = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);

	XMMATRIX world = XMMatrixScaling(1.0f, 1.0f, 1.0f);;
	world = XMMatrixMultiply(XMLoadFloat4x4(&m_constantBufferData.model), world);
	XMMATRIX view = (XMLoadFloat4x4(&m_constantBufferData.view));
	XMMATRIX proj = (XMLoadFloat4x4(&m_constantBufferData.projection));
	constbuffPerFrame.g_mViewProjection[0] = view * proj;
	constbuffPerFrame.g_mViewProjection[1] = view * proj;

	// Load shaders asynchronously.
	auto loadVSTask = DX::ReadDataAsync(L"ModelVertexShader.cso");
	auto loadPSTask = DX::ReadDataAsync(L"ModelPixelShader.cso");

	// After the vertex shader file is loaded, create the shader and input layout.
	auto createVSTask = loadVSTask.then([this](const std::vector<byte>& fileData) {
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateVertexShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_model_vertexShader
				)
			);
		static const D3D11_INPUT_ELEMENT_DESC vertexDesc[] =
		{
			{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
			{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
			{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 20, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		};

		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateInputLayout(
				vertexDesc,
				ARRAYSIZE(vertexDesc),
				&fileData[0],
				fileData.size(),
				&m_model_layout
				)
			);
	});

	// After the pixel shader file is loaded, create the shader and constant buffer.
	auto createPSTask = loadPSTask.then([this](const std::vector<byte>& fileData) {
		ID3D11Device* mDevice = m_deviceResources->GetD3DDevice();
		DX::ThrowIfFailed(
			mDevice->CreatePixelShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_model_pixelShader
				)
			);

		// vertex constant buffer 
		CD3D11_BUFFER_DESC vertexConstantBufferDesc(sizeof(VertexConstants), D3D11_BIND_CONSTANT_BUFFER);
		DX::ThrowIfFailed(
			mDevice->CreateBuffer(
				&vertexConstantBufferDesc,
				nullptr,
				&m_model_vertexCB
				)
			);

		// pixel constant buffer
		D3D11_BUFFER_DESC cbbd;
		ZeroMemory(&cbbd, sizeof(D3D11_BUFFER_DESC));

		cbbd.Usage = D3D11_USAGE_DEFAULT;
		cbbd.ByteWidth = sizeof(cbPerFrame);
		cbbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		cbbd.CPUAccessFlags = 0;
		cbbd.MiscFlags = 0;

		DX::ThrowIfFailed(
			mDevice->CreateBuffer(&cbbd, nullptr, &m_model_pixelCB)
			);

		// rasterizer state
		D3D11_RASTERIZER_DESC cmdesc;
		ZeroMemory(&cmdesc, sizeof(D3D11_RASTERIZER_DESC));
		cmdesc.FillMode = D3D11_FILL_SOLID;
		cmdesc.CullMode = D3D11_CULL_BACK;
		cmdesc.FrontCounterClockwise = true;
		mDevice->CreateRasterizerState(&cmdesc, &CCWcullMode);

		// depth stencil desc
		D3D11_DEPTH_STENCIL_DESC dssDesc;
		ZeroMemory(&dssDesc, sizeof(D3D11_DEPTH_STENCIL_DESC));
		dssDesc.DepthEnable = true;
		dssDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
		dssDesc.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
		mDevice->CreateDepthStencilState(&dssDesc, &m_model_depthSS);

		// sampler state
		D3D11_SAMPLER_DESC sampDesc;
		ZeroMemory(&sampDesc, sizeof(sampDesc));
		sampDesc.Filter = D3D11_FILTER_ANISOTROPIC;
		sampDesc.AddressU = D3D11_TEXTURE_ADDRESS_MIRROR;
		sampDesc.AddressV = D3D11_TEXTURE_ADDRESS_MIRROR;
		sampDesc.AddressW = D3D11_TEXTURE_ADDRESS_MIRROR;
		sampDesc.ComparisonFunc = D3D11_COMPARISON_NEVER;
		sampDesc.MinLOD = 0;
		sampDesc.MaxLOD = D3D11_FLOAT32_MAX;
		DX::ThrowIfFailed(
			mDevice->CreateSamplerState(&sampDesc, &m_model_samplerState)
		);

		// resource view
		DX::ThrowIfFailed(
			DirectX::CreateDDSTextureFromFile(mDevice, L"Media\\uffizi_cross.dds", NULL, &m_model_envmapSRV)
		);
	});

	(createPSTask && createVSTask).then([this]() {
		m_modelComplete = true;
	});
}

void Sample3DSceneRenderer::InitModelIndices()
{
	m_model_indexBuffer.ReleaseAndGetAddressOf();
	D3D11_BUFFER_DESC indexBufferDesc;
	ZeroMemory(&indexBufferDesc, sizeof(indexBufferDesc));
	indexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	indexBufferDesc.ByteWidth = sizeof(DWORD) * m_model_indexCount;
	indexBufferDesc.BindFlags = D3D11_BIND_INDEX_BUFFER;
	indexBufferDesc.CPUAccessFlags = 0;
	indexBufferDesc.MiscFlags = 0;

	D3D11_SUBRESOURCE_DATA iinitData;
	iinitData.pSysMem = &model_indices[0];
	m_deviceResources->GetD3DDevice()->CreateBuffer(&indexBufferDesc, &iinitData, &m_model_indexBuffer);
}

void Sample3DSceneRenderer::InitModelVertices()
{
	// m_model_vertexBuffer.ReleaseAndGetAddressOf();
	D3D11_BUFFER_DESC vertexBufferDesc;
	ZeroMemory(&vertexBufferDesc, sizeof(vertexBufferDesc));
	vertexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	vertexBufferDesc.ByteWidth = sizeof(mVertex) * model_vertices.size();
	vertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vertexBufferDesc.CPUAccessFlags = 0;
	vertexBufferDesc.MiscFlags = 0;

	D3D11_SUBRESOURCE_DATA vertexBufferData;
	ZeroMemory(&vertexBufferData, sizeof(vertexBufferData));
	vertexBufferData.pSysMem = &model_vertices[0];
	m_deviceResources->GetD3DDevice()->CreateBuffer(&vertexBufferDesc, &vertexBufferData, &m_model_vertexBuffer);
}

void Fluid::Sample3DSceneRenderer::RenderModel()
{
	if (!m_modelComplete)
	{
		return;
	}

	auto context = m_deviceResources->GetD3DDeviceContext();
	context->UpdateSubresource(m_model_pixelCB.Get(), 0, NULL, &constbuffPerFrame, 0, 0);
	context->UpdateSubresource( m_model_vertexCB.Get(), 0, NULL, &m_vertexConstants, 0, 0);

	UINT stride = sizeof(mVertex);
	UINT offset = 0;
	context->IASetVertexBuffers(
		0,
		1,
		m_model_vertexBuffer.GetAddressOf(),
		&stride,
		&offset
		);

	context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	context->IASetInputLayout(m_model_layout.Get());
	context->VSSetShader(m_model_vertexShader.Get(), nullptr, 0);
	context->VSSetConstantBuffers(0, 1, m_model_vertexCB.GetAddressOf());
	context->PSSetConstantBuffers(0, 1, m_model_pixelCB.GetAddressOf());
	context->PSSetShader(m_model_pixelShader.Get(), nullptr, 0);
	context->RSSetState(CCWcullMode.Get());
	context->OMSetDepthStencilState(m_model_depthSS.Get(), 0);
	context->PSSetShaderResources(0, 1, m_model_envmapSRV.GetAddressOf());
	//context->PSSetShaderResources(1, 1, m_model_causticSRV.GetAddressOf());
	context->Draw( model_vertices.size(), 0);
}