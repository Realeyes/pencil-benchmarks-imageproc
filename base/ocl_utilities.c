#include <stdio.h>
#include <stdlib.h>

#if defined(__APPLE__)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

/* Return the OpenCL error string for a given error number.
 */
const char *opencl_error_string(cl_int error)
{
	int errorCount;
	int index;

	static const char *errorString[] = {
		[CL_SUCCESS] = "CL_SUCCESS",
		[-CL_DEVICE_NOT_FOUND] = "CL_DEVICE_NOT_FOUND",
		[-CL_DEVICE_NOT_AVAILABLE] = "CL_DEVICE_NOT_AVAILABLE",
		[-CL_COMPILER_NOT_AVAILABLE] = "CL_COMPILER_NOT_AVAILABLE",
		[-CL_MEM_OBJECT_ALLOCATION_FAILURE] =
			"CL_MEM_OBJECT_ALLOCATION_FAILURE",
		[-CL_OUT_OF_RESOURCES] = "CL_OUT_OF_RESOURCES",
		[-CL_OUT_OF_HOST_MEMORY] = "CL_OUT_OF_HOST_MEMORY",
		[-CL_PROFILING_INFO_NOT_AVAILABLE] =
			"CL_PROFILING_INFO_NOT_AVAILABLE",
		[-CL_MEM_COPY_OVERLAP] = "CL_MEM_COPY_OVERLAP",
		[-CL_IMAGE_FORMAT_MISMATCH] = "CL_IMAGE_FORMAT_MISMATCH",
		[-CL_IMAGE_FORMAT_NOT_SUPPORTED] =
			"CL_IMAGE_FORMAT_NOT_SUPPORTED",
		[-CL_BUILD_PROGRAM_FAILURE] = "CL_BUILD_PROGRAM_FAILURE",
		[-CL_MAP_FAILURE] = "CL_MAP_FAILURE",
		[-CL_INVALID_VALUE] = "CL_INVALID_VALUE",
		[-CL_INVALID_DEVICE_TYPE] = "CL_INVALID_DEVICE_TYPE",
		[-CL_INVALID_PLATFORM] = "CL_INVALID_PLATFORM",
		[-CL_INVALID_DEVICE] = "CL_INVALID_DEVICE",
		[-CL_INVALID_CONTEXT] = "CL_INVALID_CONTEXT",
		[-CL_INVALID_QUEUE_PROPERTIES] = "CL_INVALID_QUEUE_PROPERTIES",
		[-CL_INVALID_COMMAND_QUEUE] = "CL_INVALID_COMMAND_QUEUE",
		[-CL_INVALID_HOST_PTR] = "CL_INVALID_HOST_PTR",
		[-CL_INVALID_MEM_OBJECT] = "CL_INVALID_MEM_OBJECT",
		[-CL_INVALID_IMAGE_FORMAT_DESCRIPTOR] =
			"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
		[-CL_INVALID_IMAGE_SIZE] = "CL_INVALID_IMAGE_SIZE",
		[-CL_INVALID_SAMPLER] = "CL_INVALID_SAMPLER",
		[-CL_INVALID_BINARY] = "CL_INVALID_BINARY",
		[-CL_INVALID_BUILD_OPTIONS] = "CL_INVALID_BUILD_OPTIONS",
		[-CL_INVALID_PROGRAM] = "CL_INVALID_PROGRAM",
		[-CL_INVALID_PROGRAM_EXECUTABLE] =
			"CL_INVALID_PROGRAM_EXECUTABLE",
		[-CL_INVALID_KERNEL_NAME] = "CL_INVALID_KERNEL_NAME",
		[-CL_INVALID_KERNEL_DEFINITION] =
			"CL_INVALID_KERNEL_DEFINITION",
		[-CL_INVALID_KERNEL] = "CL_INVALID_KERNEL",
		[-CL_INVALID_ARG_INDEX] = "CL_INVALID_ARG_INDEX",
		[-CL_INVALID_ARG_VALUE] = "CL_INVALID_ARG_VALUE",
		[-CL_INVALID_ARG_SIZE] = "CL_INVALID_ARG_SIZE",
		[-CL_INVALID_KERNEL_ARGS] = "CL_INVALID_KERNEL_ARGS",
		[-CL_INVALID_WORK_DIMENSION] = "CL_INVALID_WORK_DIMENSION",
		[-CL_INVALID_WORK_GROUP_SIZE] = "CL_INVALID_WORK_GROUP_SIZE",
		[-CL_INVALID_WORK_ITEM_SIZE] = "CL_INVALID_WORK_ITEM_SIZE",
		[-CL_INVALID_GLOBAL_OFFSET] = "CL_INVALID_GLOBAL_OFFSET",
		[-CL_INVALID_EVENT_WAIT_LIST] = "CL_INVALID_EVENT_WAIT_LIST",
		[-CL_INVALID_EVENT] = "CL_INVALID_EVENT",
		[-CL_INVALID_OPERATION] = "CL_INVALID_OPERATION",
		[-CL_INVALID_GL_OBJECT] = "CL_INVALID_GL_OBJECT",
		[-CL_INVALID_BUFFER_SIZE] = "CL_INVALID_BUFFER_SIZE",
		[-CL_INVALID_MIP_LEVEL] = "CL_INVALID_MIP_LEVEL",
		[-CL_INVALID_GLOBAL_WORK_SIZE] = "CL_INVALID_GLOBAL_WORK_SIZE",
		[-CL_INVALID_PROPERTY] = "CL_INVALID_PROPERTY"
	};

	errorCount = sizeof(errorString) / sizeof(errorString[0]);
	index = -error;

	return (index >= 0 && index < errorCount) ?
		errorString[index] : "Unspecified Error";
}

/* <sys/time.h> needed for time measurements.  */
#include <sys/time.h>

/* Find a GPU or a CPU associated with the first available platform.
 * If use_gpu is set, then this function first tries to look for a GPU
 * in the first available platform.
 * If this fails or if use_gpu is not set, then it tries to use the CPU.
 */
cl_device_id opencl_create_device(int use_gpu)
{
	cl_platform_id platform;
	cl_device_id dev;
	int err;

	err = clGetPlatformIDs(1, &platform, NULL);
	if (err < 0) {
		fprintf(stderr, "Error %s while looking for a platform.\n",
				opencl_error_string(err));
		exit(1);
	}

	err = CL_DEVICE_NOT_FOUND;
	if (use_gpu)
		err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &dev,
				NULL);
	if (err == CL_DEVICE_NOT_FOUND)
		err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &dev,
				NULL);
	if (err < 0) {
		fprintf(stderr, "Error %s while looking for a device.\n",
				opencl_error_string(err));
		exit(1);
	}
	return dev;
}

/* Create an OpenCL program from a source file and compile it.
 */
cl_program opencl_build_program(cl_context ctx, cl_device_id dev,
	const char* filename, const char* opencl_options)
{
	cl_program program;
	FILE *program_file;
	char *program_source, *program_log;
	size_t program_size, log_size, read;
	int err;

	program_file = fopen(filename, "r");
	if (program_file == NULL) {
		fprintf(stderr, "Could not find the source file.\n");
		exit(1);
	}
	fseek(program_file, 0, SEEK_END);
	program_size = ftell(program_file);
	rewind(program_file);
	program_source = (char *) malloc(program_size + 1);
	program_source[program_size] = '\0';
	read = fread(program_source, sizeof(char), program_size, program_file);
	if (read != program_size) {
		fprintf(stderr, "Error while reading the kernel.\n");
		exit(1);
	}
	fclose(program_file);

	program = clCreateProgramWithSource(ctx, 1,
		(const char **)&program_source, &program_size, &err);
	if (err < 0) {
		fprintf(stderr, "Could not create the program\n");
		exit(1);
	}
	free(program_source);
	err = clBuildProgram(program, 0, NULL, opencl_options, NULL, NULL);
	if (err < 0) {
		fprintf(stderr, "Could not build the program.\n");
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 0,
				NULL, &log_size);
		program_log = (char *) malloc(log_size + 1);
		program_log[log_size] = '\0';
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
				log_size + 1, program_log, NULL);
		fprintf(stderr, "%s\n", program_log);
		free(program_log);
		exit(1);
	}
	return program;
}

/* Return current time (number of seconds since the Epoch).  */
double ppcg_get_time()
{
	struct timeval Tp;
	int stat;

	stat = gettimeofday (&Tp, NULL);
	if (stat != 0)
		fprintf(stderr, "Error return from gettimeofday: %d", stat);
	return Tp.tv_sec + Tp.tv_usec * 1.0e-6;
}
