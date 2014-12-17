#ifndef HOG_PENCIL_H
#define HOG_PENCIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <pencil.h>
#include <stdint.h>
#include <stdbool.h>

#define NUMBER_OF_CELLS 1
#define NUMBER_OF_BINS 8
#define GAUSSIAN_WEIGHTS 1
#define SPARTIAL_WEIGHTS 0
#define SIGNED_HOG 1

void pencil_hog( const int rows
               , const int cols
               , const int step
               , const uint8_t image[]
               , const int num_locations
               , const float location[][2]
               , const float block_size[][2]
               , float hist[]    //out
               );

#ifdef __cplusplus
} // extern "C"
#endif

#endif //HOG_PENCIL_H