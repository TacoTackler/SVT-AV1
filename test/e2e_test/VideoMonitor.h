#ifndef _VIDEO_MONITOR_
#define _VIDEO_MONITOR_
#ifdef ENABLE_DEBUG_MONITOR
#include "SDL.h"
#include "EbSvtAv1Enc.h"

class VideoMonitor {
  public:
    VideoMonitor(const uint32_t width, const uint32_t height,
                 const uint32_t luma_stride, const uint8_t bit_depth,
                 const bool packed_ten_bit_mode, const char *name);

    ~VideoMonitor();
    void draw_frame(const uint8_t *luma, const uint8_t *cb, const uint8_t *cr);

  private:
    static uint32_t ref_cout;
    SDL_Renderer *renderer;
    SDL_Texture *texture;
    SDL_Window *window;
    uint8_t *monitor_buffer;
    const uint32_t width_;
    const uint32_t height_;
    const uint32_t luma_stride_;
    const uint8_t bit_depth_;
    const bool svt_compressed_2bit_plane;
    bool libaom_hack;
};
#endif
#endif  //_VIDEO_MONITOR_
