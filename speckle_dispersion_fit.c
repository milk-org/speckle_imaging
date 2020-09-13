/**
 * @file    speckle_dispersion_fit.c
 * @brief   Fit speckle cloud chromatic dispersion
 *
 *
 */


#include "CommandLineInterface/CLIcore.h"


// ==========================================
// Forward declaration(s)
// ==========================================

errno_t speckle_dispersion_fit(
    char *image_name
);



// ==========================================
// Command line interface wrapper function(s)
// ==========================================


static errno_t speckle_dispersion_fit__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_IMG)           
            == 0)
    {
        milk_module_example__stream_process_loop_simple(
            data.cmdargtoken[1].val.string);

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}





// ==========================================
// Register CLI command(s)
// ==========================================


errno_t speckle_dispersion_fit_addCLIcmd()
{

    RegisterCLIcommand(
        "speckledispfit",                                                                    
        __FILE__,                                                                     
        speckle_dispersion_fit__cli,
        "Fit speckle cloud chromatic dispersion",
        "<imageinput>",                  
        "speckledispfit im1",
        "milk_mspeckle_dispersion_fit(char *image_name)");

    return RETURN_SUCCESS;
}











/**
 * @brief Solve for speckle cloud chromatic dispersion
 * 
 * 
 * 
 */
errno_t speckle_dispersion_fit(
    char *image_name
)
{
	imageID IDin;
	
	IDin = image_ID(image_name);
	
	return RETURN_SUCCESS;
}
