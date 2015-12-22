TextureCube SkyMap;

struct SKYMAP_VS_OUTPUT    //output structure for skymap vertex shader
{
	float4 Pos : SV_POSITION;
	float3 texCoord : TEXCOORD;
};

SKYMAP_VS_OUTPUT SKYMAP_VS(float3 inPos : POSITION, float2 inTexCoord : TEXCOORD, float3 normal : NORMAL)
{
	SKYMAP_VS_OUTPUT output = (SKYMAP_VS_OUTPUT)0;

	output.Pos = mul(float4(inPos, 1.0f), WVP).xyww;

	output.texCoord = inPos;

	return output;        //send color and position to pixel shader
}

float4 SKYMAP_PS(SKYMAP_VS_OUTPUT input) : SV_Target
{
	return SkyMap.Sample(TriLinearSample, input.texCoord);
}

RasterizerState NoCulling
{
	CullMode = None;
};

DepthStencilState LessEqualDSS
{
	DepthFunc = LESS_EQUAL;
};

technique10 SkyMapTech
{
	pass P0
	{
		SetVertexShader(CompileShader(vs_4_0, SKYMAP_VS()));
		SetPixelShader(CompileShader(ps_4_0, SKYMAP_PS()));

		SetRasterizerState(NoCulling);
		SetDepthStencilState(LessEqualDSS, 0);
	}
}