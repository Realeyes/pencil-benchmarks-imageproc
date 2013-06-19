#include "errors.hpp"

const std::map<int, std::string> carp::opencl::errors = []()->std::map<int, std::string> {
    std::map<int, std::string> result;
    result[CL_SUCCESS] = "CL_SUCCESS";
    result[CL_DEVICE_NOT_FOUND] = "CL_DEVICE_NOT_FOUND";
    result[CL_DEVICE_NOT_AVAILABLE] = "CL_DEVICE_NOT_AVAILABLE";
    result[CL_COMPILER_NOT_AVAILABLE] = "CL_COMPILER_NOT_AVAILABLE";
    result[CL_MEM_OBJECT_ALLOCATION_FAILURE] = "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    result[CL_OUT_OF_RESOURCES] = "CL_OUT_OF_RESOURCES";
    result[CL_OUT_OF_HOST_MEMORY] = "CL_OUT_OF_HOST_MEMORY";
    result[CL_PROFILING_INFO_NOT_AVAILABLE] = "CL_PROFILING_INFO_NOT_AVAILABLE";
    result[CL_MEM_COPY_OVERLAP] = "CL_MEM_COPY_OVERLAP";
    result[CL_IMAGE_FORMAT_MISMATCH] = "CL_IMAGE_FORMAT_MISMATCH";
    result[CL_IMAGE_FORMAT_NOT_SUPPORTED] = "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    result[CL_BUILD_PROGRAM_FAILURE] = "CL_BUILD_PROGRAM_FAILURE";
    result[CL_MAP_FAILURE] = "CL_MAP_FAILURE";
    result[CL_MISALIGNED_SUB_BUFFER_OFFSET] = "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    result[CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST] = "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
#if __OPENCL_VERSION__ >= 120
    result[CL_COMPILE_PROGRAM_FAILURE] = "CL_COMPILE_PROGRAM_FAILURE";
    result[CL_LINKER_NOT_AVAILABLE] = "CL_LINKER_NOT_AVAILABLE";
    result[CL_LINK_PROGRAM_FAILURE] = "CL_LINK_PROGRAM_FAILURE";
    result[CL_DEVICE_PARTITION_FAILED] = "CL_DEVICE_PARTITION_FAILED";
    result[CL_KERNEL_ARG_INFO_NOT_AVAILABLE] = "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
    result[CL_INVALID_VALUE] = "CL_INVALID_VALUE";
#endif
    result[CL_INVALID_DEVICE_TYPE] = "CL_INVALID_DEVICE_TYPE";
    result[CL_INVALID_PLATFORM] = "CL_INVALID_PLATFORM";
    result[CL_INVALID_DEVICE] = "CL_INVALID_DEVICE";
    result[CL_INVALID_CONTEXT] = "CL_INVALID_CONTEXT";
    result[CL_INVALID_QUEUE_PROPERTIES] = "CL_INVALID_QUEUE_PROPERTIES";
    result[CL_INVALID_COMMAND_QUEUE] = "CL_INVALID_COMMAND_QUEUE";
    result[CL_INVALID_HOST_PTR] = "CL_INVALID_HOST_PTR";
    result[CL_INVALID_MEM_OBJECT] = "CL_INVALID_MEM_OBJECT";
    result[CL_INVALID_IMAGE_FORMAT_DESCRIPTOR] = "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    result[CL_INVALID_IMAGE_SIZE] = "CL_INVALID_IMAGE_SIZE";
    result[CL_INVALID_SAMPLER] = "CL_INVALID_SAMPLER";
    result[CL_INVALID_BINARY] = "CL_INVALID_BINARY";
    result[CL_INVALID_BUILD_OPTIONS] = "CL_INVALID_BUILD_OPTIONS";
    result[CL_INVALID_PROGRAM] = "CL_INVALID_PROGRAM";
    result[CL_INVALID_PROGRAM_EXECUTABLE] = "CL_INVALID_PROGRAM_EXECUTABLE";
    result[CL_INVALID_KERNEL_NAME] = "CL_INVALID_KERNEL_NAME";
    result[CL_INVALID_KERNEL_DEFINITION] = "CL_INVALID_KERNEL_DEFINITION";
    result[CL_INVALID_KERNEL] = "CL_INVALID_KERNEL";
    result[CL_INVALID_ARG_INDEX] = "CL_INVALID_ARG_INDEX";
    result[CL_INVALID_ARG_VALUE] = "CL_INVALID_ARG_VALUE";
    result[CL_INVALID_ARG_SIZE] = "CL_INVALID_ARG_SIZE";
    result[CL_INVALID_KERNEL_ARGS] = "CL_INVALID_KERNEL_ARGS";
    result[CL_INVALID_WORK_DIMENSION] = "CL_INVALID_WORK_DIMENSION";
    result[CL_INVALID_WORK_GROUP_SIZE] = "CL_INVALID_WORK_GROUP_SIZE";
    result[CL_INVALID_WORK_ITEM_SIZE] = "CL_INVALID_WORK_ITEM_SIZE";
    result[CL_INVALID_GLOBAL_OFFSET] = "CL_INVALID_GLOBAL_OFFSET";
    result[CL_INVALID_EVENT_WAIT_LIST] = "CL_INVALID_EVENT_WAIT_LIST";
    result[CL_INVALID_EVENT] = "CL_INVALID_EVENT";
    result[CL_INVALID_OPERATION] = "CL_INVALID_OPERATION";
    result[CL_INVALID_GL_OBJECT] = "CL_INVALID_GL_OBJECT";
    result[CL_INVALID_BUFFER_SIZE] = "CL_INVALID_BUFFER_SIZE";
    result[CL_INVALID_MIP_LEVEL] = "CL_INVALID_MIP_LEVEL";
    result[CL_INVALID_GLOBAL_WORK_SIZE] = "CL_INVALID_GLOBAL_WORK_SIZE";
#if __OPENCL_VERSION__ >= 120
    result[CL_INVALID_PROPERTY] = "CL_INVALID_PROPERTY";
    result[CL_INVALID_IMAGE_DESCRIPTOR] = "CL_INVALID_IMAGE_DESCRIPTOR";
    result[CL_INVALID_COMPILER_OPTIONS] = "CL_INVALID_COMPILER_OPTIONS";
    result[CL_INVALID_LINKER_OPTIONS] = "CL_INVALID_LINKER_OPTIONS";
    result[CL_INVALID_DEVICE_PARTITION_COUNT] = "CL_INVALID_DEVICE_PARTITION_COUNT";
#endif
    return result;
}(); // cl_errors
