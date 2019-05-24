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
uniform mediump float alphaLimb;
uniform mediump float betaLimb;
uniform mediump float limbGamma;
uniform bool applyToneMapping;
uniform bool scaleLimb;



void main()
{
    vec3 textureColor = vec3(texture(ourTexture, TexCoord).r);

    mediump vec3 scaledRGB = alpha * textureColor.rgb + beta;

    if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
    {
        scaledRGB.rgb = vec3(0);
    }

    vec2 pixelXY = abs(TexCoord.xy - vec2(0.5, 0.5));

    //if (radiusXYnorm != vec2(0, 0))
    //{
        //float circleEllipseEquation = pow(pixelXY.x / radiusXYnorm.x, 2) + pow(pixelXY.y / radiusXYnorm.y, 2);
//        if (circleEllipseEquation > 1)
//        {
//            if (applyToneMapping)
//            {
//                vec3 mu3 = vec3(mu);
//                /// Inverse Gaussian
//                scaledRGB = iMax * sqrt(lambda / (2.0*3.1415 * pow(scaledRGB, vec3(3.0)))) * exp(-lambda * pow(scaledRGB - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * scaledRGB)) + scaledRGB;
//                mediump vec3 gammaVec = vec3(limbGamma);
//                scaledRGB = pow(scaledRGB, gammaVec);
//            }
//            else if (scaleLimb)
//            {
//                scaledRGB = alphaLimb * textureColor.rgb + betaLimb;
//                mediump vec3 gammaVec = vec3(limbGamma);
//                scaledRGB = pow(scaledRGB, gammaVec);

//                if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
//                {
//                    scaledRGB.rgb = vec3(0);
//                }

//            }

//        }

        if (applyToneMapping)
        {
            vec3 mu3 = vec3(mu);
            /// Inverse Gaussian
            scaledRGB = iMax * sqrt(lambda / (2.0*3.1415 * pow(scaledRGB, vec3(3.0)))) * exp(-lambda * pow(scaledRGB - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * scaledRGB)) + scaledRGB;
            mediump vec3 gammaVec = vec3(limbGamma);
            scaledRGB = pow(scaledRGB, gammaVec);
        }
        else if (scaleLimb)
        {
            scaledRGB = alphaLimb * textureColor.rgb + betaLimb;
            mediump vec3 gammaVec = vec3(limbGamma);
            scaledRGB = pow(scaledRGB, gammaVec);

            if (scaledRGB.r < 0 || scaledRGB.g < 0 || scaledRGB.b < 0)
            {
                scaledRGB.rgb = vec3(0);
            }

        }
        else
        {
            mediump vec3 gammaVec = vec3(gamma);
            scaledRGB = pow(scaledRGB, gammaVec);
        }
    //}

    if (!scaleLimb)
    {
        mediump vec3 gammaVec = vec3(gamma);
        scaledRGB = pow(scaledRGB, gammaVec);
    }


//    if (applyToneMapping)
//    {
//        vec3 mu3 = vec3(mu);
//        /// Inverse Gaussian
//        scaledRGB = iMax * sqrt(lambda / (2.0*3.1415 * pow(scaledRGB, vec3(3.0)))) * exp(-lambda * pow(scaledRGB - mu3, vec3(2.0)) / (2.0 * pow(mu3, vec3(2.0)) * scaledRGB)) + scaledRGB;
//    }

    scaledRGB = scaledRGB * wbRGB;
    color = vec4(scaledRGB, 1.0);
}
