#pragma once

#include "..\Common\DeviceResources.h"
#include "ShaderStructures.h"
#include "..\Common\StepTimer.h"
#include "..\Sph\ParticleSystem.h"
#include "SurfaceGenerator.h"

namespace Fluid
{

	// This sample renderer instantiates a basic rendering pipeline.
	class Sample3DSceneRenderer
	{
	public:
		Sample3DSceneRenderer(const std::shared_ptr<DX::DeviceResources>& deviceResources);
		void CreateDeviceDependentResources();
		void CreateWindowSizeDependentResources();
		void ReleaseDeviceDependentResources();
		void Update(DX::StepTimer const& timer);
		void Render();
		void StartTracking();
		void TrackingUpdate(float positionX);
		void StopTracking();
		bool IsTracking() { return m_tracking; }

	private:
		void Rotate(float radians);
		void RenderModel();
		void LoadModel(std::wstring file_name, std::vector<DWORD> &model_indices, std::vector<mVertex> &model_vertices);
		void InitModelShaders();
		void InitModelIndices();
		void InitModelVertices();

		void RenderEnvMap();
		void RenderCaustic();
		void InitEnvMapShaders();
		void InitCausticTex();

	private:
		// Cached pointer to device resources.
		std::shared_ptr<DX::DeviceResources> m_deviceResources;

		// Direct3D resources for cube geometry.
		Microsoft::WRL::ComPtr<ID3D11InputLayout>	m_inputLayout;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_vertexBuffer;
		

		//Microsoft::WRL::ComPtr<ID3D11Buffer>		m_indexBuffer;
		Microsoft::WRL::ComPtr<ID3D11VertexShader>	m_vertexShader;
		Microsoft::WRL::ComPtr<ID3D11PixelShader>	m_pixelShader;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_constantBuffer;

		// envmap
		Microsoft::WRL::ComPtr<ID3D11VertexShader>  m_envmap_vertexShader;
		Microsoft::WRL::ComPtr<ID3D11PixelShader>   m_envmap_pixelShader;
		Microsoft::WRL::ComPtr<ID3D11InputLayout>   m_envmap_layout;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_envmap_vertexCB;  //m_pCB
		//Microsoft::WRL::ComPtr<ID3D11Buffer>		m_envmap_pixelCB;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_pVB;
		Microsoft::WRL::ComPtr<ID3D11SamplerState>	m_envmap_samplerState;
		Microsoft::WRL::ComPtr<ID3D11DepthStencilState>	m_envmap_depthSS;
		Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> m_envMapSRV;

		// model
		Microsoft::WRL::ComPtr<ID3D11InputLayout>	m_model_layout;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_model_vertexBuffer;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_model_indexBuffer;
		Microsoft::WRL::ComPtr<ID3D11VertexShader>	m_model_vertexShader;
		Microsoft::WRL::ComPtr<ID3D11PixelShader>	m_model_pixelShader;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_model_vertexCB;
		Microsoft::WRL::ComPtr<ID3D11Buffer>		m_model_pixelCB;
		Microsoft::WRL::ComPtr<ID3D11Texture2D>		m_model_caustic;
		Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> m_model_envmapSRV;
		Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> m_model_causticSRV;
		Microsoft::WRL::ComPtr<ID3D11RenderTargetView>   m_model_causticRTV;
		
		Microsoft::WRL::ComPtr<ID3D11SamplerState>	m_model_samplerState;
		Microsoft::WRL::ComPtr<ID3D11DepthStencilState> m_model_depthSS;
		Microsoft::WRL::ComPtr<ID3D11RasterizerState> CCWcullMode;
		int m_model_indexCount;

		cbPerFrame constbuffPerFrame;

		std::vector<DWORD> model_indices;
		std::vector<mVertex> model_vertices;
		std::vector<mVertex> all_vertices;

		// System resources for cube geometry.
		ModelViewProjectionConstantBuffer	m_constantBufferData;
		VertexConstants m_vertexConstants;

		uint32	m_indexCount;
        ParticleSystem m_system;
		SurfaceGenerator m_surface;
        std::vector<VertexPositionColor> m_points;
		// Variables used with the rendering loop.
		bool	m_loadingComplete;
		bool    m_modelComplete = false;
		bool	m_envmapComplete = false;
		float	m_degreesPerSecond;
		bool	m_tracking;
	};
}

