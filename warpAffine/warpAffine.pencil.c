#include "warpAffine.pencil.h"
#include <pencil.h>

static void affine( const int src_rows, const int src_cols, const int src_step, const float src[static const restrict src_rows][src_step]
                  , const int dst_rows, const int dst_cols, const int dst_step,       float dst[static const restrict dst_rows][dst_step]
                  , const float a00, const float a01, const float a10, const float a11, const float b00, const float b10
                  )
{
#pragma scop
    __pencil_assume(src_rows >  0);
    __pencil_assume(src_cols >  0);
    __pencil_assume(src_step >= src_cols);
    __pencil_assume(dst_rows >  0);
    __pencil_assume(dst_cols >  0);
    __pencil_assume(dst_step >= dst_cols);
    #pragma pencil independent
    for ( int n_r=0; n_r<dst_rows; n_r++ )
    {
        #pragma pencil independent
        for ( int n_c=0; n_c<dst_cols; n_c++ )
        {
            float o_r = a11 * n_r + a10 * n_c + b00;
            float o_c = a01 * n_r + a00 * n_c + b10;

            float r = o_r - floor(o_r);
            float c = o_c - floor(o_c);

            int coord_00_r = floor(o_r);
            int coord_00_c = floor(o_c);
            int coord_01_r = coord_00_r;
            int coord_01_c = coord_00_c + 1;
            int coord_10_r = coord_00_r + 1;
            int coord_10_c = coord_00_c;
            int coord_11_r = coord_00_r + 1;
            int coord_11_c = coord_00_c + 1;

            coord_00_r = clamp(coord_00_r, 0, src_rows);
            coord_00_c = clamp(coord_00_c, 0, src_cols);
            coord_01_r = clamp(coord_01_r, 0, src_rows);
            coord_01_c = clamp(coord_01_c, 0, src_cols);
            coord_10_r = clamp(coord_10_r, 0, src_rows);
            coord_10_c = clamp(coord_10_c, 0, src_cols);
            coord_11_r = clamp(coord_11_r, 0, src_rows);
            coord_11_c = clamp(coord_11_c, 0, src_cols);

            float A00 = src[coord_00_r][coord_00_c];
            float A10 = src[coord_10_r][coord_10_c];
            float A01 = src[coord_01_r][coord_01_c];
            float A11 = src[coord_11_r][coord_11_c];

            dst[n_r][n_c] = mix( mix(A00, A10, r), mix(A01, A11, r), c);
        }
    }
#pragma endscop
}

void pencil_affine_linear( const int src_rows, const int src_cols, const int src_step, const float src[]
                         , const int dst_rows, const int dst_cols, const int dst_step,       float dst[]
                         , const float a00, const float a01, const float a10, const float a11, const float b00, const float b10
                         )
{
    affine( src_rows, src_cols, src_step, (const float(*)[src_step])src
          , dst_rows, dst_cols, dst_step, (      float(*)[dst_step])dst
          , a00, a01, a10, a11, b00, b10
          );
}
