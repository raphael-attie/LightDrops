#version 330 core

in mediump vec2 TexCoord;

out vec4 color;

//uniform mediump sampler2D ourTexture;
uniform usampler2D ourTexture;
uniform mediump float alpha;
uniform mediump float beta;
uniform mediump float gamma;
uniform mediump vec3 wbRGB;
uniform mediump float iMax;
uniform mediump float lambda;
uniform mediump float mu;
uniform bool applyToneMapping;

void main()
{
    vec3 textureColor = vec3(texture(ourTexture, TexCoord).r);
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
        if (applyToneMapping)
        {
            /// Inverse Gaussian
            vec3 mu3 = vec3(mu);
            textureColor = iMax * sqrt(lambda / (2.0*3.1415 * pow(textureColor, vec3(3.0)))) * exp(-lambda * pow(textureColor - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * textureColor)) + textureColor;
            /// Sigmoid-like function
            //textureColor = 255.0 * textureColor / (vec3(sigmoid) + textureColor);
        }

        scaledRGB = alpha * (textureColor * wbRGB).rgb + beta;

        if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
        {
            scaledRGB.rgb = vec3(0);
        }

        mediump vec3 gammaVec = vec3(gamma);
        scaledRGB = pow(scaledRGB, gammaVec);


    }

    color = vec4(scaledRGB, 1.0);
}
