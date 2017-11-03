#version 330 core

//uniform float uKernel[9];
//uniform sampler2D uSampler;
//uniform vec2 uTextureSize;


in mediump vec2 TexCoord;

out vec4 color;

//uniform mediump sampler2D ourTexture;
uniform usampler2D ourTexture;
uniform mediump float alpha;
uniform mediump float beta;
uniform mediump float gamma;
uniform mediump vec3 wbRGB;


void main()
{
//    ivec2 texSize = textureSize(ourTexture, 0);
//    vec2 stepSize = 1.0/(uTextureSize);
//    vec3 sum = vec3(0.0);

    vec3 textureColor = texture(ourTexture, TexCoord).rgb;
    textureColor = textureColor * wbRGB;
    mediump vec3 scaledRGB = alpha * textureColor.rgb + beta;
    if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
    {
        scaledRGB.rgb = vec3(0);
    }
    mediump vec3 gammaVec = vec3(gamma);
    scaledRGB = pow(scaledRGB, gammaVec);
    color = vec4(scaledRGB, 1.0);
}
