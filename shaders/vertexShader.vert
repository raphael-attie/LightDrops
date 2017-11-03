#version 330 core

layout (location = 0) in mediump vec3 position;
layout (location = 1) in mediump vec2 texCoord;

out mediump vec2 TexCoord;

void main()
{
    gl_Position = vec4(position, 1.0f);
    TexCoord = texCoord;
}
