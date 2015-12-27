#include "pch.h"
#include "Sample3DSceneRenderer.h"
#include "..\Common\DirectXHelper.h"
#include <vector>
#include <fstream>

#include "External\CubeMap\Effects.h"
#include "External\CubeMap\Vertex.h"
#include "External\Common\MathHelper.h"
// #include <Winuser.h>

using namespace Fluid;

using namespace DirectX;
using namespace Windows::Foundation;

// Loads vertex and pixel shaders from files and instantiates the cube geometry.
Sample3DSceneRenderer::Sample3DSceneRenderer(const std::shared_ptr<DX::DeviceResources>& deviceResources) :
	m_loadingComplete(true),
	m_degreesPerSecond(45),
	m_indexCount(0),
	m_tracking(false),
	m_deviceResources(deviceResources),
    m_system(deviceResources->GetD3DDevice())
{
	// CreateDeviceDependentResources();
	CreateWindowSizeDependentResources();
	// LoadModel(L"fish.obj", model_indices, model_vertices);
	m_surface.Init(&m_system);
	// InitModelShaders();
	// init skybox
	Effects::InitAll(m_deviceResources->GetD3DDevice());
	InputLayouts::InitAll(m_deviceResources->GetD3DDevice());
	mSky = new Sky(m_deviceResources->GetD3DDevice(), L"External/CubeMap/Textures/grasscube1024.dds", 5000.0f);
	// mCam.SetPosition(0.0f, 2.0f, -2.0f);
	mCam.LookAt(XMFLOAT3(0.0f, 0.7f, -1.5f), XMFLOAT3(0.0f, -0.1f, 0.0f), XMFLOAT3(0.0f, 1.0f, 0.0f));
	// BuildSkullGeometryBuffers();
	m_system.InitParticles();
	m_points.resize(m_system.Points().extent[0]);
	std::transform(m_system.GetParticles().begin(),
		m_system.GetParticles().end(),
		m_points.begin(),
		[](const Particle& p) -> VertexPositionColor {
		return{ p.position,{ 1.0f, 1.0f, 1.0f } };
	});
}

// Initializes view parameters when the window size changes.
void Sample3DSceneRenderer::CreateWindowSizeDependentResources()
{
	Size outputSize = m_deviceResources->GetOutputSize();
	float aspectRatio = outputSize.Width / outputSize.Height;
	float fovAngleY = 70.0f * XM_PI / 180.0f;

	// This is a simple example of change that can be made when the app is in
	// portrait or snapped view.
	if (aspectRatio < 1.0f)
	{
		fovAngleY *= 2.0f;
	}

	// Note that the OrientationTransform3D matrix is post-multiplied here
	// in order to correctly orient the scene to match the display orientation.
	// This post-multiplication step is required for any draw calls that are
	// made to the swap chain render target. For draw calls to other targets,
	// this transform should not be applied.

	// This sample makes use of a right-handed coordinate system using row-major matrices.
	XMMATRIX perspectiveMatrix = XMMatrixPerspectiveFovRH(
		fovAngleY,
		aspectRatio,
		0.01f,
		100.0f
		);

	XMFLOAT4X4 orientation = m_deviceResources->GetOrientationTransform3D();

	XMMATRIX orientationMatrix = XMLoadFloat4x4(&orientation);

	XMStoreFloat4x4(
		&m_constantBufferData.projection,
		XMMatrixTranspose(perspectiveMatrix * orientationMatrix)
		);

	// Eye is at (0,0.7,1.5), looking at point (0,-0.1,0) with the up-vector along the y-axis.
	static const XMVECTORF32 eye = { 0.0f, 0.7f, 1.5f, 0.0f };
	static const XMVECTORF32 at = { 0.0f, -0.1f, 0.0f, 0.0f };
	static const XMVECTORF32 up = { 0.0f, 1.0f, 0.0f, 0.0f };

	XMStoreFloat4x4(&m_constantBufferData.view, XMMatrixTranspose(XMMatrixLookAtRH(eye, at, up)));
}

// Called once per frame, rotates the cube and calculates the model and view matrices.
void Sample3DSceneRenderer::Update(DX::StepTimer const& timer)
{
	if (!m_loadingComplete || !m_modelComplete)
		return;
	if (!m_tracking)
	{
		// Convert degrees to radians, then convert seconds to rotation angle
		float radiansPerSecond = XMConvertToRadians(m_degreesPerSecond);
		double totalRotation = timer.GetTotalSeconds() * radiansPerSecond;
		float radians = static_cast<float>(fmod(totalRotation, XM_2PI));

		Rotate(radians);
	}

    m_system.Update();

    std::vector<Concurrency::graphics::float_3> list = m_system.Points();
    std::transform(list.begin(),
        list.end(),
        m_points.begin(),
        [](const Concurrency::graphics::float_3& p) -> VertexPositionColor {
        return{ { p.x, p.y, p.z }, { 1.0f, 1.0f, 1.0f } };
    });
	// m_surface.Generate(model_vertices);
	// InitModelVertices();
	m_surface.Generate(fluid_vertices);

	int vcount = fluid_vertices.size();
	D3D11_BUFFER_DESC vbd;
	vbd.Usage = D3D11_USAGE_IMMUTABLE;
	vbd.ByteWidth = sizeof(Vertex::Basic32) * vcount;
	vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vbd.CPUAccessFlags = 0;
	vbd.MiscFlags = 0;
	D3D11_SUBRESOURCE_DATA vinitData;
	vinitData.pSysMem = &fluid_vertices[0];
	HR(m_deviceResources->GetD3DDevice()->CreateBuffer(&vbd, &vinitData, &m_fluid_vertexBuffer));

	// mCam.RotateY(0.01);
	XMMATRIX world = XMMatrixScaling(1.0f, 1.0f, 1.0f);
	world = XMMatrixMultiply(XMLoadFloat4x4(&m_constantBufferData.model), world);
	XMMATRIX view = (XMLoadFloat4x4(&m_constantBufferData.view));
	XMMATRIX proj = (XMLoadFloat4x4(&m_constantBufferData.projection));
	m_vertexConstants.MVP = XMMatrixMultiply(proj, XMMatrixMultiply(view, world));
	m_vertexConstants.World = world;
}

// Rotate the 3D cube model a set amount of radians.
void Sample3DSceneRenderer::Rotate(float radians)
{
	// Prepare to pass the updated model matrix to the shader
	XMStoreFloat4x4(&m_constantBufferData.model, XMMatrixTranspose(XMMatrixRotationY(radians)));
}

void Fluid::Sample3DSceneRenderer::RenderModel()
{
	if (!m_modelComplete)
	{
		return;
	}

	auto context = m_deviceResources->GetD3DDeviceContext();

	context->UpdateSubresource(m_model_pixelCB.Get(), 0, NULL, &constbuffPerFrame, 0, 0);
	
	// Prepare the constant buffer to send it to the graphics device.
	context->UpdateSubresource(
		m_model_vertexCB.Get(),
		0,
		NULL,
		&m_vertexConstants,
		0,
		0
		);

	// Each vertex is one instance of the VertexPositionColor struct.
	UINT stride = sizeof(mVertex);
	UINT offset = 0;

	//UINT stride = sizeof(Concurrency::graphics::float_3);
	//UINT offset = 0;
	context->IASetVertexBuffers(
		0,
		1,
		m_model_vertexBuffer.GetAddressOf(),
		&stride,
		&offset
		);
	/*context->IASetIndexBuffer(
		m_model_indexBuffer.Get(),
		DXGI_FORMAT_R32_UINT, // Each index is one 16-bit unsigned integer (short).
		0
		);*/

	context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	context->IASetInputLayout(m_model_layout.Get());

	// Attach our vertex shader.
	context->VSSetShader(
		m_model_vertexShader.Get(),
		nullptr,
		0
		);

	// Send the constant buffer to the graphics device.
	context->VSSetConstantBuffers(
		0,
		1,
		m_model_vertexCB.GetAddressOf()
		);

	context->PSSetConstantBuffers(0, 1, m_model_pixelCB.GetAddressOf());

	// Attach our pixel shader.
	context->PSSetShader(
		m_model_pixelShader.Get(),
		nullptr,
		0
		);
	context->RSSetState(CCWcullMode.Get());
	context->OMSetDepthStencilState(DSLessEqual.Get(), 0);

	// Draw the objects.
	context->Draw(
		model_vertices.size(),
		0
		);
}

void Fluid::Sample3DSceneRenderer::LoadModel(std::wstring file_name, std::vector<DWORD> &model_indices, std::vector<mVertex> &model_vertices)
{
	std::wifstream fileIn(file_name);

	std::vector<XMFLOAT3> vertPos;
	std::vector<XMFLOAT3> vertNorm;
	std::vector<XMFLOAT2> vertTexCoord;

	std::vector<int> vertPosIndex;
	std::vector<int> vertNormIndex;
	std::vector<int> vertTCIndex;

	wchar_t checkChar;
	while (fileIn.good()) {
		checkChar = fileIn.get();
		switch (checkChar) {
		case 'v':
			checkChar = fileIn.get();
			if (checkChar == ' ') {
				float vx, vy, vz;
				fileIn >> vx >> vy >> vz;
				vertPos.push_back(XMFLOAT3(vx, vy, vz));
			}
			else if (checkChar == 't') {
				float vtu, vtv;
				fileIn >> vtu >> vtv;
				vertTexCoord.push_back(XMFLOAT2(vtu, vtv));
			}
			else if (checkChar == 'n') {
				float vnx, vny, vnz;
				fileIn >> vnx >> vny >> vnz;
				vertNorm.push_back(XMFLOAT3(vnx, vny, vnz));
			}
			break;
		case 'f':
			for (int i = 0; i < 3; i++) {
				unsigned int vi, ti, ni; wchar_t c;
				fileIn >> vi >> c >> ti >> c >> ni;
				int index_count = vertPos.size();
				int j = 0, n = vertPosIndex.size();
				for (; j < n; j++) {
					if (vertPosIndex[j] == vi && vertNormIndex[j] == ni && vertTCIndex[j] == ti) {
						model_indices.push_back(j);
						break;
					}
				}
				if (j >= n) {
					vertPosIndex.push_back(vi); vertNormIndex.push_back(ni); vertTCIndex.push_back(ti);
					mVertex mv;
					mv.pos = vertPos[vi-1]; mv.normal = vertNorm[ni-1]; mv.tex = vertTexCoord[ti-1];
					model_vertices.push_back(mv);
					model_indices.push_back(n);
				}
				//mVertex mv;
				//mv.pos = vertPos[vi - 1]; mv.normal = vertNorm[ni - 1]; mv.tex = vertTexCoord[ti - 1];
				//all_vertices.push_back(mv);
			}
		default:
			while (checkChar != '\n' && fileIn.good()) checkChar = fileIn.get();
			break;
		}
	}

	m_model_indexCount = model_indices.size();

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
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreatePixelShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_model_pixelShader
				)
			);

		CD3D11_BUFFER_DESC vertexConstantBufferDesc(sizeof(VertexConstants), D3D11_BIND_CONSTANT_BUFFER);
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateBuffer(
				&vertexConstantBufferDesc,
				nullptr,
				&m_model_vertexCB
				)
			);

		D3D11_BUFFER_DESC cbbd;
		ZeroMemory(&cbbd, sizeof(D3D11_BUFFER_DESC));

		cbbd.Usage = D3D11_USAGE_DEFAULT;
		cbbd.ByteWidth = sizeof(cbPerFrame);
		cbbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		cbbd.CPUAccessFlags = 0;
		cbbd.MiscFlags = 0;

		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateBuffer(&cbbd, nullptr, &m_model_pixelCB)
			);

		D3D11_RASTERIZER_DESC cmdesc;
		ZeroMemory(&cmdesc, sizeof(D3D11_RASTERIZER_DESC));
		cmdesc.FillMode = D3D11_FILL_SOLID;
		cmdesc.CullMode = D3D11_CULL_BACK;
		cmdesc.FrontCounterClockwise = true;
		m_deviceResources->GetD3DDevice()->CreateRasterizerState(&cmdesc, &CCWcullMode);

		D3D11_DEPTH_STENCIL_DESC dssDesc;
		ZeroMemory(&dssDesc, sizeof(D3D11_DEPTH_STENCIL_DESC));
		dssDesc.DepthEnable = true;
		dssDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
		dssDesc.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;

		m_deviceResources->GetD3DDevice()->CreateDepthStencilState(&dssDesc, &DSLessEqual);
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

void Sample3DSceneRenderer::StartTracking()
{
	m_tracking = true;

}

// When tracking, the 3D cube can be rotated around its Y axis by tracking pointer position relative to the output screen width.
void Sample3DSceneRenderer::TrackingUpdate(float positionX)
{
	if (m_tracking)
	{
		float radians = XM_2PI * 2.0f * positionX / m_deviceResources->GetOutputSize().Width;
		Rotate(radians);
	}
}

void Sample3DSceneRenderer::StopTracking()
{
	m_tracking = false;
}

// Renders one frame using the vertex and pixel shaders.
void Sample3DSceneRenderer::Render()
{
	// Loading is asynchronous. Only draw geometry after it's loaded.
	if (!m_loadingComplete)
	{
		return;
	}

	auto context = m_deviceResources->GetD3DDeviceContext();

	
	mCam.UpdateViewMatrix();
	mSky->Draw(context, mCam);
	/*
	// Prepare the constant buffer to send it to the graphics device.
	context->UpdateSubresource(
		m_constantBuffer.Get(),
		0,
		NULL,
		&m_constantBufferData,
		0,
		0
		);

	// Each vertex is one instance of the VertexPositionColor struct.
	UINT stride = sizeof(VertexPositionColor);
    UINT offset = 0;
    D3D11_MAPPED_SUBRESOURCE mappedResource;
    DX::ThrowIfFailed(
        context->Map(
            m_vertexBuffer.Get(), 0, D3D11_MAP_WRITE_NO_OVERWRITE, 0, &mappedResource)
        );
    memcpy(mappedResource.pData, m_points.data(), m_points.size() * sizeof(m_points[0]));
    context->Unmap(m_vertexBuffer.Get(), 0);

    //UINT stride = sizeof(Concurrency::graphics::float_3);
    //UINT offset = 0;
    context->IASetVertexBuffers(
        0,
        1,
        m_vertexBuffer.GetAddressOf(),
        &stride,
        &offset
        );
	//context->IASetIndexBuffer(
	//	m_indexBuffer.Get(),
	//	DXGI_FORMAT_R16_UINT, // Each index is one 16-bit unsigned integer (short).
	//	0
	//	);

	context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);

	context->IASetInputLayout(m_inputLayout.Get());

	// Attach our vertex shader.
	context->VSSetShader(
		m_vertexShader.Get(),
		nullptr,
		0
		);

	// Send the constant buffer to the graphics device.
	context->VSSetConstantBuffers(
		0,
		1,
		m_constantBuffer.GetAddressOf()
		);

	// Attach our pixel shader.
	context->PSSetShader(
		m_pixelShader.Get(),
		nullptr,
		0
		);

	// Draw the objects.
	context->Draw(
		m_points.size(),
		0
		);*/

	RenderFluid();
	
}

void Sample3DSceneRenderer::CreateDeviceDependentResources()
{
	// Load shaders asynchronously.
	auto loadVSTask = DX::ReadDataAsync(L"SampleVertexShader.cso");
	auto loadPSTask = DX::ReadDataAsync(L"SamplePixelShader.cso");

	// After the vertex shader file is loaded, create the shader and input layout.
	auto createVSTask = loadVSTask.then([this](const std::vector<byte>& fileData) {
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateVertexShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_vertexShader
				)
			);
        m_system.InitParticles();
        m_points.resize(m_system.Points().extent[0]);
        std::transform(m_system.GetParticles().begin(),
            m_system.GetParticles().end(),
            m_points.begin(),
            [](const Particle& p) -> VertexPositionColor {
            return { p.position, { 1.0f, 1.0f, 1.0f } };
        });
		static const D3D11_INPUT_ELEMENT_DESC vertexDesc [] =
		{
			{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
			{ "COLOR", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		};

		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateInputLayout(
				vertexDesc,
				ARRAYSIZE(vertexDesc),
				&fileData[0],
				fileData.size(),
				&m_inputLayout
				)
			);
	});

	// After the pixel shader file is loaded, create the shader and constant buffer.
	auto createPSTask = loadPSTask.then([this](const std::vector<byte>& fileData) {
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreatePixelShader(
				&fileData[0],
				fileData.size(),
				nullptr,
				&m_pixelShader
				)
			);

		CD3D11_BUFFER_DESC constantBufferDesc(sizeof(ModelViewProjectionConstantBuffer) , D3D11_BIND_CONSTANT_BUFFER);
		DX::ThrowIfFailed(
			m_deviceResources->GetD3DDevice()->CreateBuffer(
				&constantBufferDesc,
				nullptr,
				&m_constantBuffer
				)
			);
	});

	// Once both shaders are loaded, create the mesh.
	auto createCubeTask = (createPSTask && createVSTask).then([this] () {

		// Load mesh vertices. Each vertex has a position and a color.
		//static const VertexPositionColor cubeVertices[] = 
		//{
		//	{XMFLOAT3(-0.5f, -0.5f, -0.5f), XMFLOAT3(0.0f, 0.0f, 0.0f)},
		//	{XMFLOAT3(-0.5f, -0.5f,  0.5f), XMFLOAT3(0.0f, 0.0f, 1.0f)},
		//	{XMFLOAT3(-0.5f,  0.5f, -0.5f), XMFLOAT3(0.0f, 1.0f, 0.0f)},
		//	{XMFLOAT3(-0.5f,  0.5f,  0.5f), XMFLOAT3(0.0f, 1.0f, 1.0f)},
		//	{XMFLOAT3( 0.5f, -0.5f, -0.5f), XMFLOAT3(1.0f, 0.0f, 0.0f)},
		//	{XMFLOAT3( 0.5f, -0.5f,  0.5f), XMFLOAT3(1.0f, 0.0f, 1.0f)},
		//	{XMFLOAT3( 0.5f,  0.5f, -0.5f), XMFLOAT3(1.0f, 1.0f, 0.0f)},
		//	{XMFLOAT3( 0.5f,  0.5f,  0.5f), XMFLOAT3(1.0f, 1.0f, 1.0f)},
		//};

		//D3D11_SUBRESOURCE_DATA vertexBufferData = {0};
		//vertexBufferData.pSysMem = cubeVertices;
		//vertexBufferData.SysMemPitch = 0;
		//vertexBufferData.SysMemSlicePitch = 0;
		//CD3D11_BUFFER_DESC vertexBufferDesc(sizeof(cubeVertices), D3D11_BIND_VERTEX_BUFFER);
		//DX::ThrowIfFailed(
		//	m_deviceResources->GetD3DDevice()->CreateBuffer(
		//		&vertexBufferDesc,
		//		&vertexBufferData,
		//		&m_vertexBuffer
		//		)
		//	);

        //DX::ThrowIfFailed(
        //    m_system.GetDataBuffer(reinterpret_cast<void**>(m_vertexBuffer.GetAddressOf()))
        //    );
        D3D11_SUBRESOURCE_DATA vertexBufferData = { 0 };
        vertexBufferData.pSysMem = m_points.data();
        vertexBufferData.SysMemPitch = 0;
        vertexBufferData.SysMemSlicePitch = 0;
        CD3D11_BUFFER_DESC vertexBufferDesc;
        vertexBufferDesc.Usage = D3D11_USAGE_DYNAMIC;
        vertexBufferDesc.ByteWidth = sizeof(VertexPositionColor) * m_points.size();
        vertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
        vertexBufferDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
        vertexBufferDesc.MiscFlags = 0;
        vertexBufferDesc.StructureByteStride = 0;
        DX::ThrowIfFailed(
            m_deviceResources->GetD3DDevice()->CreateBuffer(
                &vertexBufferDesc,
                &vertexBufferData,
                &m_vertexBuffer)
            );

		// Load mesh indices. Each trio of indices represents
		// a triangle to be rendered on the screen.
		// For example: 0,2,1 means that the vertices with indexes
		// 0, 2 and 1 from the vertex buffer compose the 
		// first triangle of this mesh.
		static const unsigned short cubeIndices [] =
		{
			0,2,1, // -x
			1,2,3,

			4,5,6, // +x
			5,7,6,

			0,1,5, // -y
			0,5,4,

			2,6,7, // +y
			2,7,3,

			0,4,6, // -z
			0,6,2,

			1,3,7, // +z
			1,7,5,
		};

		m_indexCount = ARRAYSIZE(cubeIndices);

		//D3D11_SUBRESOURCE_DATA indexBufferData = {0};
		//indexBufferData.pSysMem = cubeIndices;
		//indexBufferData.SysMemPitch = 0;
		//indexBufferData.SysMemSlicePitch = 0;
		//CD3D11_BUFFER_DESC indexBufferDesc(sizeof(cubeIndices), D3D11_BIND_INDEX_BUFFER);
		//DX::ThrowIfFailed(
		//	m_deviceResources->GetD3DDevice()->CreateBuffer(
		//		&indexBufferDesc,
		//		&indexBufferData,
		//		&m_indexBuffer
		//		)
		//	);
	});

	// Once the cube is loaded, the object is ready to be rendered.
	createCubeTask.then([this] () {
		m_loadingComplete = true;
	});
}

void Sample3DSceneRenderer::ReleaseDeviceDependentResources()
{
	m_loadingComplete = false;
	m_vertexShader.Reset();
	m_inputLayout.Reset();
	m_pixelShader.Reset();
	m_constantBuffer.Reset();
	m_vertexBuffer.Reset();
	//m_indexBuffer.Reset();
}

/*
void Sample3DSceneRenderer::BuildSkullGeometryBuffers()
{
	std::ifstream fin("skull.txt");

	if (!fin)
	{
		// MessageBox(0, L"Models/skull.txt not found.", 0, 0);
		return;
	}

	UINT vcount = 0;
	UINT tcount = 0;
	std::string ignore;

	fin >> ignore >> vcount;
	fin >> ignore >> tcount;
	fin >> ignore >> ignore >> ignore >> ignore;

	std::vector<Vertex::Basic32> vertices(vcount);
	for (UINT i = 0; i < vcount; ++i)
	{
		fin >> vertices[i].Pos.x >> vertices[i].Pos.y >> vertices[i].Pos.z;
		fin >> vertices[i].Normal.x >> vertices[i].Normal.y >> vertices[i].Normal.z;
	}

	fin >> ignore;
	fin >> ignore;
	fin >> ignore;

	mSkullIndexCount = 3 * tcount;
	std::vector<UINT> indices(mSkullIndexCount);
	for (UINT i = 0; i < tcount; ++i)
	{
		fin >> indices[i * 3 + 0] >> indices[i * 3 + 1] >> indices[i * 3 + 2];
	}

	fin.close();

	D3D11_BUFFER_DESC vbd;
	vbd.Usage = D3D11_USAGE_IMMUTABLE;
	vbd.ByteWidth = sizeof(Vertex::Basic32) * vcount;
	vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vbd.CPUAccessFlags = 0;
	vbd.MiscFlags = 0;
	D3D11_SUBRESOURCE_DATA vinitData;
	vinitData.pSysMem = &vertices[0];
	HR(m_deviceResources->GetD3DDevice()->CreateBuffer(&vbd, &vinitData, &mSkullVB));

	//
	// Pack the indices of all the meshes into one index buffer.
	//

	D3D11_BUFFER_DESC ibd;
	ibd.Usage = D3D11_USAGE_IMMUTABLE;
	ibd.ByteWidth = sizeof(UINT) * mSkullIndexCount;
	ibd.BindFlags = D3D11_BIND_INDEX_BUFFER;
	ibd.CPUAccessFlags = 0;
	ibd.MiscFlags = 0;
	D3D11_SUBRESOURCE_DATA iinitData;
	iinitData.pSysMem = &indices[0];
	HR(m_deviceResources->GetD3DDevice()->CreateBuffer(&ibd, &iinitData, &mSkullIB));
}
*/

void Sample3DSceneRenderer::RenderFluid()
{
	auto md3dImmediateContext = m_deviceResources->GetD3DDeviceContext();
	// md3dImmediateContext->ClearDepthStencilView(mDepthStencilView, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);

	md3dImmediateContext->IASetInputLayout(InputLayouts::Basic32);
	md3dImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	UINT stride = sizeof(Vertex::Basic32);
	UINT offset = 0;

	mCam.UpdateViewMatrix();

	XMMATRIX view = mCam.View();
	XMMATRIX proj = mCam.Proj();
	XMMATRIX viewProj = mCam.ViewProj();

	float blendFactor[] = { 0.0f, 0.0f, 0.0f, 0.0f };

	// Set per frame constants.
	Effects::BasicFX->SetDirLights(mDirLights);
	Effects::BasicFX->SetEyePosW(mCam.GetPosition());
	Effects::BasicFX->SetCubeMap(mSky->CubeMapSRV());

	// Figure out which technique to use.  Skull does not have texture coordinates,
	// so we need a separate technique for it, and not every surface is reflective,
	// so don't pay for cubemap look up.

	ID3DX11EffectTechnique* activeTexTech = Effects::BasicFX->Light1TexTech;
	ID3DX11EffectTechnique* activeReflectTech = Effects::BasicFX->Light1TexReflectTech;
	ID3DX11EffectTechnique* activeSkullTech = Effects::BasicFX->Light1ReflectTech;
	switch (mLightCount)
	{
	case 1:
		activeTexTech = Effects::BasicFX->Light1TexTech;
		activeReflectTech = Effects::BasicFX->Light1TexReflectTech;
		activeSkullTech = Effects::BasicFX->Light1ReflectTech;
		break;
	case 2:
		activeTexTech = Effects::BasicFX->Light2TexTech;
		activeReflectTech = Effects::BasicFX->Light2TexReflectTech;
		activeSkullTech = Effects::BasicFX->Light2ReflectTech;
		break;
	case 3:
		activeTexTech = Effects::BasicFX->Light3TexTech;
		activeReflectTech = Effects::BasicFX->Light3TexReflectTech;
		activeSkullTech = Effects::BasicFX->Light3ReflectTech;
		break;
	}

	XMMATRIX world;
	XMMATRIX worldInvTranspose;
	XMMATRIX worldViewProj;

	D3DX11_TECHNIQUE_DESC techDesc;
	activeSkullTech->GetDesc(&techDesc);
	for (UINT p = 0; p < techDesc.Passes; ++p)
	{
		// Draw the skull.

		md3dImmediateContext->IASetVertexBuffers(0, 1, &m_fluid_vertexBuffer, &stride, &offset);
		// md3dImmediateContext->IASetIndexBuffer(mSkullIB, DXGI_FORMAT_R32_UINT, 0);

		//XMMATRIX skullScale = XMMatrixScaling(1.0f, 1.0f, 1.0f);
		//XMMATRIX skullOffset = XMMatrixTranslation(0.0f, 1.0f, 0.0f);
		//XMFLOAT4X4 mSkullWorld;
		//XMStoreFloat4x4(&mSkullWorld, XMMatrixMultiply(skullScale, skullOffset));
		//world = XMLoadFloat4x4(&mSkullWorld);

		//XMMATRIX world = XMMatrixScaling(1.0f, 1.0f, 1.0f);
		//world = XMMatrixMultiply(XMLoadFloat4x4(&m_constantBufferData.model), world);
		//XMMATRIX view = (XMLoadFloat4x4(&m_constantBufferData.view));
		//XMMATRIX proj = (XMLoadFloat4x4(&m_constantBufferData.projection));
		//m_vertexConstants.MVP = XMMatrixMultiply(proj, XMMatrixMultiply(view, world));
		//m_vertexConstants.World = world;

		XMMATRIX world = XMMatrixScaling(0.5f, 0.5f, 0.5f);
		// world = XMMatrixMultiply(XMLoadFloat4x4(&m_constantBufferData.model), world);
		
		worldInvTranspose = MathHelper::InverseTranspose(world);
		worldViewProj = world*view*proj;

		Material mSkullMat;
		mSkullMat.Ambient = XMFLOAT4(0.2f, 0.2f, 0.2f, 1.0f);
		mSkullMat.Diffuse = XMFLOAT4(0.2f, 0.2f, 0.2f, 1.0f);
		mSkullMat.Specular = XMFLOAT4(0.8f, 0.8f, 0.8f, 16.0f);
		mSkullMat.Reflect = XMFLOAT4(0.5f, 0.5f, 0.5f, 1.0f);

		Effects::BasicFX->SetWorld(world);
		Effects::BasicFX->SetWorldInvTranspose(worldInvTranspose);
		Effects::BasicFX->SetWorldViewProj(worldViewProj);
		Effects::BasicFX->SetMaterial(mSkullMat);

		activeSkullTech->GetPassByIndex(p)->Apply(0, md3dImmediateContext);
		md3dImmediateContext->Draw(fluid_vertices.size(), 0);
	}

	// mSky->Draw(md3dImmediateContext, mCam);

	// restore default states, as the SkyFX changes them in the effect file.
	md3dImmediateContext->RSSetState(0);
	md3dImmediateContext->OMSetDepthStencilState(0, 0);
}

void Sample3DSceneRenderer::InitFluidVertices()
{
	int vcount = fluid_vertices.size();
	D3D11_BUFFER_DESC vbd;
	vbd.Usage = D3D11_USAGE_IMMUTABLE;
	vbd.ByteWidth = sizeof(Vertex::Basic32) * vcount;
	vbd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vbd.CPUAccessFlags = 0;
	vbd.MiscFlags = 0;
	D3D11_SUBRESOURCE_DATA vinitData;
	vinitData.pSysMem = &fluid_vertices[0];
	HR(m_deviceResources->GetD3DDevice()->CreateBuffer(&vbd, &vinitData, &m_fluid_vertexBuffer));
}