#version 330 core

//uniform float uKernel[9];
//uniform sampler2D uSampler;
//uniform vec2 uTextureSize;


in mediump vec2 TexCoord;

//out vec4 color;
out vec4 textureColor;

//uniform mediump sampler2D ourTexture;
uniform usampler2D ourTexture;
uniform mediump float alpha;
uniform mediump float beta;
uniform mediump float gamma;
uniform mediump vec3 wbRGB;

// Kernel convolution reference:
// https://www.taylorpetrick.com/blog/post/convolution-part1
// http://coding-experiments.blogspot.com/2010/07/convolution.html
void main()
{
    // Initialize kernel. Discussion on the web says this may not be compatible
    // with some openGL ES versions (e.g on ios), in which case, need to do it the old-fashioned way.
    // (see references above)
//    float[9] uKernel = float[] (-1.,-1.,-1.,
//                                -1., 9., -1.,
//                                -1., -1., -1.);

//    vec4 sum = vec4(0.0);
//    ivec2 texSize = textureSize(ourTexture, 0);
//    vec2 stepSize = 1.0/vec2(float(texSize.x), float(texSize.y));

//    sum += texture(ourTexture, vec2(TexCoord.x - stepSize.x, TexCoord.y - stepSize.y)) * uKernel[0];
//    sum += texture(ourTexture, vec2(TexCoord.x, TexCoord.y - stepSize.y)) * uKernel[1];
//    sum += texture(ourTexture, vec2(TexCoord.x + stepSize.x, TexCoord.y - stepSize.y)) * uKernel[2];
//    sum += texture(ourTexture, vec2(TexCoord.x - stepSize.x, TexCoord.y)) * uKernel[3];
//    sum += texture(ourTexture, vec2(TexCoord.x, TexCoord.y)) * uKernel[4];
//    sum += texture(ourTexture, vec2(TexCoord.x + stepSize.x, TexCoord.y)) * uKernel[5];
//    sum += texture(ourTexture, vec2(TexCoord.x - stepSize.x, TexCoord.y + stepSize.y)) * uKernel[6];
//    sum += texture(ourTexture, vec2(TexCoord.x, TexCoord.y + stepSize.y)) * uKernel[7];
//    sum += texture(ourTexture, vec2(TexCoord.x + stepSize.x, TexCoord.y + stepSize.y)) * uKernel[8];

    textureColor = texture(ourTexture, TexCoord);
    textureColor.rgb = textureColor.rgb * wbRGB;
    textureColor.rgb = alpha * textureColor.rgb + beta;
    //mediump vec3 scaledRGB = alpha * textureColor.rgb + beta;
    if (textureColor.r < 0 || textureColor.g < 0 || textureColor.b < 0)
    {
        textureColor.rgb = vec3(0);
    }
    mediump vec3 gammaVec = vec3(gamma);
    textureColor.rgb = pow(textureColor.rgb, gammaVec);
    textureColor.a = 1.0;
    //color = vec4(textureColor.rgb, 1.0);
}
