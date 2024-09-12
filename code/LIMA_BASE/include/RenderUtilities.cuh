#pragma once

#include <cuda_runtime.h>



namespace RenderUtilities {
	enum ATOM_TYPE { NONE, O, C, P, N, H, SOL, S, LIMA_CUSTOM };

    __device__ __host__ inline ATOM_TYPE RAS_getTypeFromAtomletter(char atom) {
        switch (atom)
        {
        case 'C':
            return ATOM_TYPE::C;
        case 'O':
            return ATOM_TYPE::O;
        case 'N':
            return ATOM_TYPE::N;
        case 'H':
            return ATOM_TYPE::H;
        case 'P':
            return ATOM_TYPE::P;
        case 'S':
            return ATOM_TYPE::S;
        case 'l':
            return ATOM_TYPE::LIMA_CUSTOM;
        default:
            printf("Unknown atom type: %c\n", atom);
            return ATOM_TYPE::NONE;
        }
    }

    __device__ __host__ inline float4 getColor(ATOM_TYPE atom_type) {
        switch (atom_type)
        {
        case ATOM_TYPE::SOL:
            return float4{ 0x03 / 255.0f, 0xa9 / 255.0f, 0xf4 / 255.0f, 1.0f };
        case ATOM_TYPE::H:
            return float4{ 0xF1 / 255.0f, 0xF1 / 255.0f, 0xF1 / 255.0f, 1.0f };
        case ATOM_TYPE::O:
            return float4{ 0xE0 / 255.0f, 0x20 / 255.0f, 0x20 / 255.0f, 1.0f };
        case ATOM_TYPE::C:
            return float4{ 0x20 / 255.0f, 0x10 / 255.0f, 0x90 / 255.0f, 1.0f };
        case ATOM_TYPE::P:
            return float4{ 0xFC / 255.0f, 0xF7 / 255.0f, 0x5E / 255.0f, 1.0f };
        case ATOM_TYPE::N:
            return float4{ 0x2E / 255.0f, 0x8B / 255.0f, 0x57 / 255.0f, 1.0f };
        case ATOM_TYPE::S:
            return float4{ 0xF4 / 255.0f, 0xC4 / 255.0f, 0x30 / 255.0f, 1.0f };
        case ATOM_TYPE::NONE:
            return float4{ 0xFF / 255.0f, 0x00 / 255.0f, 0xFF / 255.0f, 1.0f };
        default:
            return float4{ 0.0f, 0.0f, 0.0f, 1.0f };
        }
    }

    //https://en.wikipedia.org/wiki/Van_der_Waals_radius
    __device__ float inline getRadius(ATOM_TYPE atom_type) {
        switch (atom_type)
        {
        case ATOM_TYPE::H:
            return 0.109f * 0.25f;   // Make smaller for visibility
        case ATOM_TYPE::C:
            return 0.17f;
        case ATOM_TYPE::N:
            return 0.155f;
        case ATOM_TYPE::O:
            return 0.152f;
        case ATOM_TYPE::SOL:
            return 0.152f * 0.4f;   // Make smaller for visibility
        case ATOM_TYPE::P:
            return 0.18f;
        case ATOM_TYPE::S:
            return 0.189f;
        case ATOM_TYPE::LIMA_CUSTOM:
            return 0.1f;
        case ATOM_TYPE::NONE:
            return .5f;
        default:
            return 1.f;
        }
    }

    float4 inline HSLtoRGB(float h, float s, float b) {

        float i = std::floor(h * 6.0f);
        float f = h * 6.0f - i;
        float p = b * (1.0f - s);
        float q = b * (1.0f - f * s);
        float t = b * (1.0f - (1.0f - f) * s);

        switch (static_cast<int>(i) % 6) {
        case 0: return { b, t, p, 1.0f };
        case 1: return { q, b, p, 1.0f };
        case 2: return { p, b, t, 1.0f };
        case 3: return { p, q, b, 1.0f };
        case 4: return { t, p, b, 1.0f };
        case 5: return { b, p, q, 1.0f };
        default: return { b, p, q, 1.0f };
        }
    }

    float4 static GetColorInGradient(float value) {
        value = std::fmaxf(value, 0.f);
        value = std::fminf(value, 1.f);


        return float4(value, 0.f, 0.1f, 1.f);
        //value = std::max(value, 0.f);

        // Set constant saturation and brightness
        float saturation = 1.0f;
        float brightness = .8f;

        // Use value as hue (0-1 maps to 0-360 degrees in HSB)
        return HSLtoRGB(value, saturation, brightness);
    }


}