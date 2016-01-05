cbuffer cbPerObject : register(b0)
{
	row_major matrix   g_mWorldViewProjection	: packoffset(c0);
}

struct VS_Input
{
	float4 Pos : POSITION;
};

struct VS_Output
{
	float4 Pos : SV_POSITION;
	float3 Tex : TEXCOORD0;
};

VS_Output main( VS_Input Input )
{
	VS_Output Output;

	Output.Pos = Input.Pos;
	Output.Tex = normalize(mul(Input.Pos, g_mWorldViewProjection));

	return Output;
}