/**\file */
#ifndef SLIC_DECLARATIONS_S_LQCD_H
#define SLIC_DECLARATIONS_S_LQCD_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define S_LQCD_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_cfactor Interface Parameter "cfactor".
 * \param [in] param_ka Interface Parameter "ka".
 * \param [in] instream_times1kernel_gauge0 The stream should be of size 2985984 bytes.
 * \param [in] instream_times1kernel_gauge1 The stream should be of size 1440000 bytes.
 * \param [in] instream_times1kernel_spinor_in The stream should be of size 995328 bytes.
 * \param [out] outstream_times1kernel_spinor_out The stream should be of size 480000 bytes.
 */
void S_LQCD(
	float param_cfactor,
	float param_ka,
	const float *instream_times1kernel_gauge0,
	const float *instream_times1kernel_gauge1,
	const float *instream_times1kernel_spinor_in,
	float *outstream_times1kernel_spinor_out);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_cfactor Interface Parameter "cfactor".
 * \param [in] param_ka Interface Parameter "ka".
 * \param [in] instream_times1kernel_gauge0 The stream should be of size 2985984 bytes.
 * \param [in] instream_times1kernel_gauge1 The stream should be of size 1440000 bytes.
 * \param [in] instream_times1kernel_spinor_in The stream should be of size 995328 bytes.
 * \param [out] outstream_times1kernel_spinor_out The stream should be of size 480000 bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *S_LQCD_nonblock(
	float param_cfactor,
	float param_ka,
	const float *instream_times1kernel_gauge0,
	const float *instream_times1kernel_gauge1,
	const float *instream_times1kernel_spinor_in,
	float *outstream_times1kernel_spinor_out);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	float param_cfactor; /**<  [in] Interface Parameter "cfactor". */
	float param_ka; /**<  [in] Interface Parameter "ka". */
	const float *instream_times1kernel_gauge0; /**<  [in] The stream should be of size 2985984 bytes. */
	const float *instream_times1kernel_gauge1; /**<  [in] The stream should be of size 1440000 bytes. */
	const float *instream_times1kernel_spinor_in; /**<  [in] The stream should be of size 995328 bytes. */
	float *outstream_times1kernel_spinor_out; /**<  [out] The stream should be of size 480000 bytes. */
} S_LQCD_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void S_LQCD_run(
	max_engine_t *engine,
	S_LQCD_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_nonblock(
	max_engine_t *engine,
	S_LQCD_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void S_LQCD_run_group(max_group_t *group, S_LQCD_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_group_nonblock(max_group_t *group, S_LQCD_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void S_LQCD_run_array(max_engarray_t *engarray, S_LQCD_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *S_LQCD_run_array_nonblock(max_engarray_t *engarray, S_LQCD_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* S_LQCD_convert(max_file_t *maxfile, S_LQCD_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* S_LQCD_init(void);

/* Error handling functions */
int S_LQCD_has_errors(void);
const char* S_LQCD_get_errors(void);
void S_LQCD_clear_errors(void);
/* Free statically allocated maxfile data */
void S_LQCD_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int S_LQCD_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int S_LQCD_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_S_LQCD_H */

