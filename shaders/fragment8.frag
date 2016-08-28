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
uniform mediump vec2 radiusXYnorm;
uniform bool applyToneMapping;
uniform bool useInverseGaussian;

void main()
{
    vec3 textureColor = vec3(texture(ourTexture, TexCoord).r);

    vec2 pixelXY = abs(TexCoord.xy - vec2(0.5, 0.5));
    if (radiusXYnorm != vec2(0, 0))
    {
        float circleEllipseEquation = pow(pixelXY.x / radiusXYnorm.x, 2) + pow(pixelXY.y / radiusXYnorm.y, 2);
        if (circleEllipseEquation > 1)
        {
            if (applyToneMapping)
            {
                vec3 mu3 = vec3(mu);

                /// Inverse Gaussian
                textureColor = iMax * sqrt(lambda / (2.0*3.1415 * pow(textureColor, vec3(3.0)))) * exp(-lambda * pow(textureColor - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * textureColor)) + textureColor;
                /// Sigmoid-like function
                //textureColor = 255.0 * textureColor / (vec3(sigmoid) + textureColor);
            }
        }
    }
    else
    {
        if (applyToneMapping)
        {
            vec3 mu3 = vec3(mu);

            /// Inverse Gaussian
            textureColor = iMax * sqrt(lambda / (2.0*3.1415 * pow(textureColor, vec3(3.0)))) * exp(-lambda * pow(textureColor - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * textureColor)) + textureColor;
            /// Sigmoid-like function
            //textureColor = 255.0 * textureColor / (vec3(sigmoid) + textureColor);
        }
    }

    mediump vec3 scaledRGB;
    scaledRGB = alpha * textureColor.rgb + beta;
    scaledRGB = scaledRGB * wbRGB;

    if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
    {
        scaledRGB.rgb = vec3(0);
    }

    mediump vec3 gammaVec = vec3(gamma);
    scaledRGB = pow(scaledRGB, gammaVec);

    color = vec4(scaledRGB, 1.0);
}
