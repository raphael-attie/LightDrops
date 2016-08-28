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
    mediump vec3 scaledRGB;

    scaledRGB = (alpha * textureColor.rgb + beta)*wbRGB;
    //scaledRGB = alpha * (textureColor * wbRGB).rgb + beta;

    /// Clip minimum values to 0. Is that even necessary?
    if (scaledRGB.r < 0){ scaledRGB.r = 0; }
    if (scaledRGB.g < 0){ scaledRGB.g = 0; }
    if (scaledRGB.b < 0){ scaledRGB.b = 0; }

    mediump vec3 gammaVec = vec3(gamma);
    scaledRGB = pow(scaledRGB, gammaVec);

    color = vec4(scaledRGB, 1.0);
}
