#version 330 core

in mediump vec2 TexCoord;

out vec4 color;

//uniform mediump sampler2D ourTexture;
uniform sampler2D ourTexture;
uniform mediump float alpha;
uniform mediump float beta;
uniform mediump float gamma;
uniform mediump vec3 wbRGB;


void main()
{
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

