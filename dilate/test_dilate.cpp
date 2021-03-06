#include "utility.hpp"
#include "dilate.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <prl.h>
#include <chrono>

void time_dilate( const std::vector<carp::record_t>& pool, const std::vector<int>& elemsizes, int iteration )
{
    bool first_execution_opencv = true, first_execution_pencil = true;

    carp::Timing timing("dilate image");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & item : pool ) {
            for ( auto & elemsize : elemsizes ) {
                // acquiring the image for the test
                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Point anchor( elemsize/2, elemsize/2 );
                cv::Size ksize(elemsize, elemsize);
                cv::Mat structuring_element = cv::getStructuringElement( cv::MORPH_ELLIPSE, ksize, anchor );

                cv::Mat cpu_result, gpu_result, pen_result;
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy;

                {
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cv::dilate( cpu_gray, cpu_result, structuring_element, anchor, 1, cv::BORDER_CONSTANT );
                    const auto cpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_cpu = cpu_end - cpu_start;
                }
                {
                    // Execute the kernel at least once before starting to take time measurements so that the OpenCV kernel gets compiled. The following run is not included in time measurements.
                    if (first_execution_opencv)
                    {
                        cv::ocl::oclMat gpu_gray(cpu_gray);
                        cv::ocl::oclMat result;
                        cv::ocl::dilate( gpu_gray, result, structuring_element, anchor, 1, cv::BORDER_CONSTANT );
                        first_execution_opencv = false;
                    }

                    const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                    cv::ocl::oclMat gpu_gray(cpu_gray);
                    cv::ocl::oclMat result;
                    cv::ocl::dilate( gpu_gray, result, structuring_element, anchor, 1, cv::BORDER_CONSTANT );
                    gpu_result = result;
                    const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                    elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                }
                {
                    pen_result = cv::Mat(cpu_gray.size(), CV_8U);

                    if (first_execution_pencil)
                    {
                        pencil_dilate( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr(), pen_result.step1(), pen_result.ptr(), structuring_element.rows, structuring_element.cols, structuring_element.step1(), structuring_element.ptr(), anchor.x, anchor.y);
                        first_execution_pencil = false;
                    }

                    prl_timings_reset();
                    prl_timings_start();
                    pencil_dilate( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr()
                                 , pen_result.step1(), pen_result.ptr()
                                 , structuring_element.rows, structuring_element.cols, structuring_element.step1(), structuring_element.ptr()
                                 , anchor.x, anchor.y
                                 );
                    prl_timings_stop();
                    // Dump execution times for PENCIL code.
                    prl_timings_dump();
                }
                // Verifying the results
                if ( (cv::norm(cpu_result - gpu_result) > 0.01) || (cv::norm(cpu_result - pen_result) > 0.01) ) {
                    std::cerr << "ERROR: Results don't match. Writing calculated images." << std::endl;
                    std::cerr << "CPU norm:" << cv::norm(cpu_result) << std::endl;
                    std::cerr << "GPU norm:" << cv::norm(gpu_result) << std::endl;
                    std::cerr << "PEN norm:" << cv::norm(pen_result) << std::endl;
                    std::cerr << "GPU-CPU norm:" << cv::norm(gpu_result, cpu_result) << std::endl;
                    std::cerr << "PEN-CPU norm:" << cv::norm(pen_result, cpu_result) << std::endl;

                    cv::imwrite( "dilate_cpu.png", cpu_result );
                    cv::imwrite( "dilate_gpu.png", gpu_result );
                    cv::imwrite( "dilate_pen.png", pen_result );
                    cv::imwrite( "dilate_cpugpu.png", cv::abs(cpu_result-gpu_result) );
                    cv::imwrite( "dilate_cpupen.png", cv::abs(cpu_result-pen_result) );
                    throw std::runtime_error("The OpenCL or PENCIL results are not equivalent with the C++ results.");
                }
                // Dump execution times for OpenCV calls.
                timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy );
            }
        }
    }
}

int main(int argc, char* argv[])
{
    prl_init((prl_init_flags)(PRL_TARGET_DEVICE_DYNAMIC | PRL_PROFILING_ENABLED));

    try {
	std::cout << "This executable is iterating over all the files passed to it as an argument. " << std::endl;

	auto pool = carp::get_pool(argc, argv);

#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_dilate( pool, { 5 }, 1 );
#else
        time_dilate( pool, { 3, 5 }, 25 );
#endif

        prl_shutdown();
        return EXIT_SUCCESS;
    }catch(const std::exception& e) {
        std::cout << e.what() << std::endl;

        prl_shutdown();
        return EXIT_FAILURE;
    }
}
