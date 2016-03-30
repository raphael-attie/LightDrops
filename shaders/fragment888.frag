#version 330 core

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
    vec3 textureColor = texture(ourTexture, TexCoord).rgb;
    //textureColor = textureColor * wbRGB;
    mediump vec3 scaledRGB;

    if (textureColor.r == 255 && textureColor.g != 255)
    {
        scaledRGB = vec3(255, 0, 0);
    }
    else if (textureColor.r != 255 && textureColor.g == 255)
    {
        scaledRGB = vec3(0, 255, 0);
    }
    else
    {
        scaledRGB = alpha * (textureColor * wbRGB).rgb + beta;
        if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
        {
            scaledRGB.rgb = vec3(0);
        }
        mediump vec3 gammaVec = vec3(gamma);
        scaledRGB = pow(scaledRGB, gammaVec);
    }

    //mediump vec3 scaledRGB = textureColor.rgb;

    color = vec4(scaledRGB, 1.0);
}
