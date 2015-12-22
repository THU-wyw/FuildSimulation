

cbuffer VertexConstants
{
	float4x4 WVP;
	float4x4 World;
};

struct VS_INPUT
{
	float4 inPos : POSITION;
	float2 inTexCoord : TEXCOORD;
	float3 normal : NORMAL;
};

struct VS_OUTPUT
{
	float4 Pos : SV_POSITION;
	float4 worldPos : POSITION;
	float2 TexCoord : TEXCOORD;
	float3 normal : NORMAL;
};

VS_OUTPUT main(VS_INPUT input)
{
	VS_OUTPUT output;

	output.Pos = mul(input.inPos, WVP);
	output.worldPos = mul(input.inPos, World);

	output.normal = mul(input.normal, World);

	output.TexCoord = input.inTexCoord;

	return output;
}