
#include "VideoMonitor.h"
#include "stdio.h"
#include "string.h"
#ifdef ENABLE_DEBUG_MONITOR
uint32_t VideoMonitor::ref_cout = 0;

VideoMonitor::VideoMonitor(const uint32_t width, const uint32_t height,
                           const uint32_t luma_stride, const uint8_t bit_depth,
                           const bool use_compressed_2bit_plane,
                           const char *name)
    : width_(width),
      height_(height),
      luma_stride_(luma_stride),
      bit_depth_(bit_depth),
      svt_compressed_2bit_plane(use_compressed_2bit_plane) {
#ifdef ENABLE_DEBUG_MONITOR
    printf("init_monitor:%dx%d.%d\r\n", width, height, bit_depth);
    if (ref_cout == 0) {
        if (SDL_Init(SDL_INIT_VIDEO) != 0) {
            printf("SDL_Init Error!\r\n");
            return;
        }
    }
    ref_cout++;
    int x = 0;
    int y = 20;
    libaom_hack = false;
    if (strcmp(name, "Ref decode") == 0) {
        x = 600;
        libaom_hack = true;
    }
    if (strcmp(name, "Recon") == 0) {
        x = 300;
    }

    window =
        SDL_CreateWindow(name, x, y, width, height, SDL_WINDOW_ALLOW_HIGHDPI);
    if (window == nullptr) {
        printf("SDL_CreateWindow Error!\r\n");
        return;
    }

    renderer = SDL_CreateRenderer(window, -1, 0);
    if (renderer == nullptr) {
        printf("SDL_CreateRenderer Error!\r\n");
        return;
    }
    // TODO: Check more format
    texture = SDL_CreateTexture(renderer,
                                SDL_PIXELFORMAT_IYUV,
                                SDL_TEXTUREACCESS_STREAMING,
                                width,
                                height);
    if (texture == nullptr) {
        printf("SDL_CreateTexture Error!\r\n");
        return;
    }
    monitor_buffer = (uint8_t *)SDL_calloc(1, width * height * 3);

#ifndef _WIN32
    SDL_Event e;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                quit = true;
            }
            if (e.type == SDL_KEYDOWN) {
                quit = true;
            }
            if (e.type == SDL_MOUSEBUTTONDOWN) {
                quit = true;
            }
        }
    }
#endif
#endif
}

VideoMonitor::~VideoMonitor() {
#ifdef ENABLE_DEBUG_MONITOR
    SDL_DestroyWindow(window);
    ref_cout--;
    if (ref_cout == 0)
        SDL_Quit();
#endif
}
void VideoMonitor::draw_frame(const uint8_t *luma, const uint8_t *cb,
                              const uint8_t *cr) {
#ifdef ENABLE_DEBUG_MONITOR
    const unsigned int luma_len = width_ * height_;
    if (bit_depth_ == 8 || ((bit_depth_ == 10) && svt_compressed_2bit_plane)) {
        unsigned int i = 0;
        for (uint32_t l = 0; l < height_; l++) {
            const uint8_t *p = luma + l * luma_stride_;
            for (uint32_t r = 0; r < width_; r++) {
                if (libaom_hack)
                    monitor_buffer[i++] = p[r * 2];
                else
                    monitor_buffer[i++] = p[r];
            }
        }

        for (uint32_t l = 0; l < (height_ >> 1); l++) {
            const uint8_t *p = cb + l * (luma_stride_ >> 1);
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                if (libaom_hack)
                    monitor_buffer[i++] = p[r * 2];
                else
                    monitor_buffer[i++] = p[r];
            }
        }

        for (uint32_t l = 0; l < (height_ >> 1); l++) {
            const uint8_t *p = cr + l * (luma_stride_ >> 1);
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                if (libaom_hack)
                    monitor_buffer[i++] = p[r * 2];
                else
                    monitor_buffer[i++] = p[r];
            }
        }
    } else if (bit_depth_ > 8 && !svt_compressed_2bit_plane) {
        uint32_t i = 0;
        const uint8_t bit_shift = (bit_depth_ - 8);

        for (uint32_t l = 0; l < height_; l++) {
            const uint16_t *p = (uint16_t *)(luma + l * luma_stride_);
            for (uint32_t r = 0; r < width_; r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }

        for (int l = 0; l < (height_ >> 1); l++) {
            const uint16_t *p = (uint16_t *)(cb + l * (luma_stride_ >> 1));
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }

        for (uint32_t l = 0; l < (height_ >> 1); l++) {
            const uint16_t *p = (uint16_t *)(cr + l * (luma_stride_ >> 1));
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }
    } else if (bit_depth_ > 8 && svt_compressed_2bit_plane) {
        uint32_t i = 0;
        const uint8_t bit_shift = (bit_depth_ - 8);

        for (uint32_t l = 0; l < height_; l++) {
            const uint16_t *p = (uint16_t *)(luma + l * luma_stride_);
            for (uint32_t r = 0; r < width_; r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }

        for (uint32_t l = 0; l < (height_ >> 1); l++) {
            const uint16_t *p = (uint16_t *)(cb + l * (luma_stride_ >> 1));
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }

        for (uint32_t l = 0; l < (height_ >> 1); l++) {
            const uint16_t *p = (uint16_t *)(cr + l * (luma_stride_ >> 1));
            for (uint32_t r = 0; r < (width_ >> 1); r++) {
                monitor_buffer[i++] = p[r] >> bit_shift & 0xFF;
            }
        }
    }
    SDL_UpdateTexture(texture, NULL, monitor_buffer, width_);

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);
// usleep(1000*1000);
#endif
}
#endif
