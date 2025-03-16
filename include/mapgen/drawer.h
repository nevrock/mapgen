#ifndef DRAWER_H
#define DRAWER_H

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

class Drawer {
public:
    // Uses stb_image_write_png to write a PNG image
    static void write(const char* filename, int width, int height, int channels, const unsigned char* data) {
        stbi_write_png(filename, width, height, channels, data, width * channels);
    }
};

#endif // DRAWER_H