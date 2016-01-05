TextureCube<float4> EnvMap: register(t0);
Texture2D<float4> LightMap : register(t1);
//Texture2D<float4> txReflection: register(t1);
//Texture2D<float4> txRefraction: register(t2);
//Texture2D<float4> txCaustic: register(t3);

SamplerState g_samLinear : register(s0);

struct Light
{
	float3 pos;
	float  range;
	float3 dir;
	float cone;
	float3 att;
	float4 ambient;
	float4 diffuse;
};

cbuffer cbPerFrame
{
	float4x4 g_mViewProjection[2];
	Light light;
};

struct VS_OUTPUT
{
	float4 Pos : SV_POSITION;
	float4 worldPos : POSITION;
	float4 wvPos : POSITION1;
	float2 TexCoord : TEXCOORD;
	float3 normal : NORMAL;
};

#define IOR 1.333333
#define R0Constant (((1.0- (1.0/IOR) )*(1.0- (1.0/IOR) ))/((1.0+ (1.0/IOR) )*(1.0+ (1.0/IOR) )))
#define R0Inv (1.0 - R0Constant)

float FresnelApprox(float3 incident, float3 normal)
{
	return R0Constant + R0Inv * pow(1.0 - dot(incident, normal), 5.0);
}

float3 intersectP(float3 v, float3 r)
{
	const static float3 np = { 0, 1, 0 };
	float3 vrp = (r - dot(r, np) * np) * dot(v, np);
	return float3(vrp.x + v.x, 0, vrp.z + v.z);
}

float2 intersectTex(float3 pos)
{
	float4 vpos1 = mul(float4(pos, 1.0f), g_mViewProjection[1]);
	float2 spos = vpos1.xy / vpos1.w;
	spos.y = -spos.y;
	spos = spos * 0.5f + 0.5f;
	return spos;
}

float3 SampleWorld(float3 v, float3 r)
{
	float dln = dot(r, float3(0, 1, 0));
	float c = 0;

	c = EnvMap.Sample(g_samLinear, r).xyz;
	return c;
}

const static float3 bodycolor = { 1, 1, 1 };
const static float4 eyePos = {0.0f, 0.7f, 1.5f, 0.0f};
const static float4 PlaneColor = { 0.3f, 0.5f, 1.0f, 1 };

float4 main(VS_OUTPUT input) : SV_TARGET
{
	//int2 iTex = input.Pos.xy;
	//float3 p_pos = PositionMap[iTex].xyz;
	//clip(eyePos.w - 0.5f);

	float3 norm = normalize(input.normal);
	float3 lightToPixelVec = light.pos - input.worldPos.xyz;
	float3 p_pos = input.worldPos.xyz;

	float3 eyedir = normalize(p_pos - eyePos.xyz);
	float3 refl = normalize(reflect(eyedir, norm));
	float3 refra = normalize(refract(eyedir, norm, 1 / 1.3333f));
	float wf = FresnelApprox(refl, norm);

	float3 cReflect = SampleWorld(p_pos, refl) * wf;
	float3 cRefract = SampleWorld(p_pos, refra) * (1.0f - wf);
	//float3 cReflect = { 0.3f, 0.5f, 1.0f};
	//float3 cRefract = { 0.0f, 0.0f, 0.1f};

	float3 envc = 0;
	float sp = saturate(dot(refl, lightToPixelVec));
	//sp = pow(sp, 400.0f) * 4.0f;
	sp = pow(sp, 50.0f) * 0.8f;
	
	envc += sp;

	float3 fc = envc + cReflect + cRefract * bodycolor;
	return float4(fc, 1.0f);

	//input.normal = normalize(input.normal);

	/*
	//Set diffuse color of material
	float4 diffuse = float4(0.8f, 0.0f, 0.0f, 0.8f);

	float3 finalColor = float3(0.0f, 0.0f, 0.0f);

	//Create the vector between light position and pixels position
	float3 lightToPixelVec = light.pos - input.worldPos;

	//Find the distance between the light pos and pixel pos
	float d = length(lightToPixelVec);

	//Add the ambient light
	float3 finalAmbient = diffuse * light.ambient;

	//If pixel is too far, return pixel color with ambient light
	if (d > light.range)
	return float4(finalAmbient, diffuse.a);

	//Turn lightToPixelVec into a unit length vector describing
	//the pixels direction from the lights position
	lightToPixelVec /= d;

	//Calculate how much light the pixel gets by the angle
	//in which the light strikes the pixels surface
	float howMuchLight = dot(lightToPixelVec, input.normal);

	//If light is striking the front side of the pixel
	if (howMuchLight > 0.0f)
	{
		//Add light to the finalColor of the pixel
		finalColor += diffuse * light.diffuse;

		//Calculate Light's Distance Falloff factor
		finalColor /= (light.att[0] + (light.att[1] * d)) + (light.att[2] * (d*d));

		//Calculate falloff from center to edge of pointlight cone
		finalColor *= pow(max(dot(-lightToPixelVec, light.dir), 0.0f), light.cone);
	}

	//make sure the values are between 1 and 0, and add the ambient
	finalColor = saturate(finalColor + finalAmbient);

	//Return Final Color
	return float4(finalColor, diffuse.a);
	*/
}


/*
float3 SampleWorld(float3 v, float3 r)
{
	float dln = dot(r, float3(0, 1, 0));
	float c = 0;

	if (dln < 0)
	{
		float2 cTex = intersectTex(intersectP(v, r));
		float4 cCaustic = LightMap.Sample(g_samLinear, cTex);
		if (cCaustic.w < 0.001f || dot(cCaustic.xyz, float3(1, 1, 1)) < 0.05f)
		{
			c = EnvMap.Sample(g_samLinear, r).xyz;
		}
		else
		{
			c = cCaustic.xyz;
		}
	}
	else
	{
		c = EnvMap.Sample(g_samLinear, r).xyz;
	}
	return c;
}

float2 intersectTex(float3 pos)
{
	float4 vpos1 = mul(float4(pos, 1.0f), g_mViewProjection[1]);
	float2 spos = vpos1.xy / vpos1.w;
	spos.y = -spos.y;
	spos = spos * 0.5f + 0.5f;
	return spos;
}

float3 intersectP(float3 v, float3 r)
{
	const static float3 np = { 0, 1, 0 };
	float3 vrp = (r - dot(r, np) * np) * dot(v, np);
	return float3(vrp.x + v.x, 0, vrp.z + v.z);
}

// reference: http://www.digitalrune.com/Blog/Post/1758/Water-Rendering
float ComputeCaustics(float3 position, float3 lightIntensity)
{
	// We could make the caustics depth-dependent, but without a high sample count,
	// this introduces artifacts.
	float depth = 1;  //SurfaceLevel - position.y

	float3 refractionToSurface = -DirectionalLightDirection * depth;
	float3 causticPosition = position + refractionToSurface;

	// For method 1 (Add), the intial value must be 0.
	// For method 2 (Multiply), the initial value must be 1.
	float caustics = 1;

	// Combine NxN samples.
	for (int i = 0; i < CausticsSampleCount; i++)
	{
		for (int j = 0; j < CausticsSampleCount; j++)
		{
			// A horizontal offset for the sample position
			float3 offset = depth * float3(
				(i - CausticsSampleCount / 2) * CausticsSampleOffset,
				0,
				(j - CausticsSampleCount / 2) * CausticsSampleOffset);

			// Lookup position.
			float2 texCoord = (causticPosition.xz + offset.xz) * WaveMapScale + WaveMapOffset;

			// The direction to the lookup position is:
			float3 refractedDirection = refractionToSurface + offset;

			// Sample wave normal map.
			float3 n = tex2Dlod(WaveNormalSampler, float4(texCoord, 0, 0)).xzy * 2 - 1;

			// Use wave normal to perturb refractedDirection.
			refractedDirection -= depth * CausticsDistortion * float3(n.x, 0, n.z);

			// Light contribution of this sample.
			float sampleIntensity = pow(
				max(0.0001, dot(normalize(refractedDirection), -DirectionalLightDirection)),
				CausticsExponent);

			// Method 1: Add intensity of samples. Creates smoother caustics.
			//caustics += sampleIntensity;

			// Method 2: Multiply intensities. Creates caustics with more contrast.
			caustics *= sampleIntensity;
		}
	}

	caustics *= CausticsIntensity;

	// Optional: Apply a non-linear curve to enhance contrast.
	//caustics = pow(caustics, 3);

	caustics *= lightIntensity;

	return caustics;
}
*/