/**
 * @file    speckle_dispersion_fit.c
 * @brief   Fit speckle cloud chromatic dispersion
 *
 *
 */


#include <math.h>

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"


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
        speckle_dispersion_fit(
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
    uint32_t xsize, ysize, zsize;

    FILE *fpout;


    IDin = image_ID(image_name);
    xsize = data.image[IDin].md[0].size[0];
    ysize = data.image[IDin].md[0].size[1];
    zsize = data.image[IDin].md[0].size[2];

    uint64_t xysize = xsize * ysize;


    int FASTMODE = 0; // set to 1 to skip tests and intermediate results

    // radius for PSF photocenter estimation
    // pixels further out will be ignored
    float centradius = 50.0;

    // exclude pixels inside this radius for radiation center
    float innerradius = 20.0;

    uint32_t pixgap = 1;
    float pixfluxlim = 100.0;



    uint32_t *arraysize;
    arraysize = (uint32_t *) malloc(sizeof(float) * 3);
    arraysize[0] = xsize;
    arraysize[1] = ysize;
    arraysize[2] = 14;
    imageID IDout = create_image_ID("testout", 3, arraysize, _DATATYPE_FLOAT, 0, 0);
    arraysize[2] = zsize;
    imageID IDoutCent = create_image_ID("radcent", 3, arraysize, _DATATYPE_FLOAT, 0,
                                        0);
    free(arraysize);



    // active pixels in image (bright pixels)
    uint32_t *iiarray;
    uint32_t *jjarray;
    uint64_t *iijjarray;
    float *pvarray;
    long NBactpix = 0;


    // solution
    float flux = 0.0;
    float photcentx = 0.0;
    float photcenty = 0.0;
    float radcentx = 0.0;
    float radcenty = 0.0;



    fpout = fopen("radcent.dat", "w");


    for(uint32_t kk = 0; kk < zsize; kk++)
    {
        printf("Processing image %4u / %4u\n", kk, zsize);

        double xcent0 = 0.0;
        double ycent0 = 0.0;
        double fcent0 = 0.0;


        // zero
        if(FASTMODE == 0)
        {
            for(uint32_t ii = 0; ii < xysize; ii++)
            {
                data.image[IDout].array.F[9 * xysize + ii] = 0.0;
                data.image[IDout].array.F[10 * xysize + ii] = 0.0;
                data.image[IDout].array.F[11 * xysize + ii] = 0.0;
                data.image[IDout].array.F[12 * xysize + ii] = 0.0;
            }
        }


        // find out coarse PSF center
        NBactpix = 0;
        for(uint32_t ii = 0; ii < xsize; ii++)
        {
            for(uint32_t jj = 0; jj < ysize; jj++)
            {
                float pixval = data.image[IDin].array.F[kk * xysize + jj * xsize + ii];
                xcent0 += pixval * ii;
                ycent0 += pixval * jj;
                fcent0 += pixval;

                if(pixval > pixfluxlim)
                {
                    NBactpix ++;
                }
            }
        }
        xcent0 /= fcent0;
        ycent0 /= fcent0;


        iiarray = (uint32_t *) malloc(sizeof(uint32_t) * NBactpix);
        jjarray = (uint32_t *) malloc(sizeof(uint32_t) * NBactpix);
        iijjarray = (uint64_t *) malloc(sizeof(uint64_t) * NBactpix);
        pvarray = (float *) malloc(sizeof(float) * NBactpix);


        double xcent1 = 0.0;
        double ycent1 = 0.0;
        double fcent1 = 0.0;

        uint32_t iimin, iimax, jjmin, jjmax;
        iimin = pixgap;
        iimax = xsize - pixgap;
        jjmin = pixgap;
        jjmax = ysize - pixgap;

        if(xcent0 - centradius > pixgap)
        {
            iimin = (uint32_t)(xcent0 - centradius);
        }
        if(xcent0 + centradius < xsize - pixgap)
        {
            iimax = (uint32_t)(xcent0 + centradius);
        }

        if(ycent0 - centradius > pixgap)
        {
            jjmin = (uint32_t)(ycent0 - centradius);
        }
        if(ycent0 + centradius < xsize - pixgap)
        {
            jjmax = (uint32_t)(ycent0 + centradius);
        }


        long actpix = 0;
        for(uint32_t ii = iimin; ii < iimax; ii++)
        {
            for(uint32_t jj = jjmin; jj < jjmax; jj++)
            {
                float pixval = data.image[IDin].array.F[kk * xysize + jj * xsize + ii];

                xcent1 += pixval * ii;
                ycent1 += pixval * jj;
                fcent1 += pixval;


                if(pixval > pixfluxlim)
                {
                    iiarray[actpix] = ii;
                    jjarray[actpix] = jj;
                    iijjarray[actpix] = jj * xsize + ii;
                    pvarray[actpix] = pixval;
                    actpix++;

                    if(FASTMODE == 0)
                    {

                        data.image[IDout].array.F[jj * xsize + ii] = pixval;

                        float pixvalmm = data.image[IDin].array.F[kk * xysize + (jj - 1) * xsize +
                                         (ii - 1)];
                        float pixvalm0 = data.image[IDin].array.F[kk * xysize + (jj - 1) * xsize +
                                         (ii)];
                        float pixvalmp = data.image[IDin].array.F[kk * xysize + (jj - 1) * xsize +
                                         (ii + 1)];
                        float pixval0m = data.image[IDin].array.F[kk * xysize + (jj) * xsize +
                                         (ii - 1)];
                        float pixval0p = data.image[IDin].array.F[kk * xysize + (jj) * xsize +
                                         (ii + 1)];
                        float pixvalpm = data.image[IDin].array.F[kk * xysize + (jj + 1) * xsize +
                                         (ii - 1)];
                        float pixvalp0 = data.image[IDin].array.F[kk * xysize + (jj + 1) * xsize +
                                         (ii)];
                        float pixvalpp = data.image[IDin].array.F[kk * xysize + (jj + 1) * xsize +
                                         (ii + 1)];

                        float pixvalsum = 1.0 * (pixval + pixvalmm + pixvalm0 + pixvalmp + pixval0m +
                                                 pixval0p + pixvalpm + pixvalp0 + pixvalpp);

                        data.image[IDout].array.F[1 * xysize + jj * xsize + ii] =
                            (pixval - pixvalmm) / pixvalsum;
                        data.image[IDout].array.F[2 * xysize + jj * xsize + ii] =
                            (pixval - pixvalm0) / pixvalsum;
                        data.image[IDout].array.F[3 * xysize + jj * xsize + ii] =
                            (pixval - pixvalmp) / pixvalsum;
                        data.image[IDout].array.F[4 * xysize + jj * xsize + ii] =
                            (pixval - pixval0m) / pixvalsum;
                        data.image[IDout].array.F[5 * xysize + jj * xsize + ii] =
                            (pixval - pixval0p) / pixvalsum;
                        data.image[IDout].array.F[6 * xysize + jj * xsize + ii] =
                            (pixval - pixvalpm) / pixvalsum;
                        data.image[IDout].array.F[7 * xysize + jj * xsize + ii] =
                            (pixval - pixvalp0) / pixvalsum;
                        data.image[IDout].array.F[8 * xysize + jj * xsize + ii] =
                            (pixval - pixvalpp) / pixvalsum;
                    }

                }
            }
        }
        NBactpix = actpix;
        xcent1 /= fcent1;
        ycent1 /= fcent1;


        printf("    Photocenter X = %12.3f   -> %12.3f\n", xcent0, xcent1);
        printf("    Photocenter Y = %12.3f   -> %12.3f\n", ycent0, ycent1);
        printf("    Total Flux    = %12.4g   -> %12.4g\n", fcent0, fcent1);

        printf("    %6lu pixels processed\n", actpix);


        photcentx = xcent1;
        photcenty = ycent1;



        float radcentval = 1.0e30; // initialization
        radcentx = 0.0;
        radcenty = 0.0;


        uint32_t iircmin, iircmax, jjrcmin, jjrcmax;
        uint32_t rcstep = 4;

        float rcradius = centradius * 0.7;

        iircmin = (int)(xcent1 - rcradius);
        iircmax = (int)(xcent1 + rcradius);

        jjrcmin = (int)(ycent1 - rcradius);
        jjrcmax = (int)(ycent1 + rcradius);

        if(iircmin < 0)
        {
            iircmin = 0;
        }
        if(iircmax > xsize)
        {
            iircmax = xsize;
        }

        if(jjrcmin < 0)
        {
            jjrcmin = 0;
        }

        if(jjrcmax > ysize)
        {
            jjrcmax = ysize;
        }




        float radfact = 1.1;


        int OKloop = 1;
        while(OKloop == 1)
        {

            for(uint32_t iirc = iircmin; iirc < iircmax; iirc += rcstep)
            {
                for(uint32_t jjrc = jjrcmin; jjrc < jjrcmax; jjrc += rcstep)
                {

                    float dxr = 1.0 * iirc - xcent1;
                    float dyr = 1.0 * jjrc - ycent1;

                    if((dxr * dxr + dyr * dyr) < rcradius * rcradius)
                    {
                        double errfit1 = 0.0;
                        double errfit2 = 0.0;
                        long errfitcnt = 0;

                        for(actpix = 0; actpix < NBactpix; actpix++)
                        {
                            uint32_t ii = iiarray[actpix];
                            uint32_t jj = jjarray[actpix];
                            uint64_t iijj = iijjarray[actpix];
                            float pixval = pvarray[actpix];

                            // radius from photocenter
                            float dxp = 1.0 * ii - xcent1;
                            float dyp = 1.0 * jj - ycent1;
                            float pradius = sqrt(dxp * dxp + dyp * dyp);

                            if(pradius > innerradius)
                            {


                                int dii = ii - iirc;
                                int djj = jj - jjrc;

                                // offset to center pixel
                                float dx = 1.0 * dii;
                                float dy = 1.0 * djj;

                                // offset for radially offset pixel
                                float dx1 = dx * radfact;
                                float dy1 = dy * radfact;

                                // pixel motion
                                float ddx1 = dx1 - dx;
                                float ddy1 = dy1 - dy;

                                //
                                float dx2 = dx - ddy1;
                                float dy2 = dy + ddx1;

                                float ii1f = 1.0 * iirc + dx1;
                                float jj1f = 1.0 * jjrc + dy1;

                                float ii2f = 1.0 * iirc + dx2;
                                float jj2f = 1.0 * jjrc + dy2;


                                int ii1 = (int)(ii1f);
                                int jj1 = (int)(jj1f);

                                int ii2 = (int)(ii2f);
                                int jj2 = (int)(jj2f);


                                float ii1frac = ii1f - ii1;
                                float jj1frac = jj1f - jj1;

                                float ii2frac = ii2f - ii2;
                                float jj2frac = jj2f - jj2;


                                if((ddx1 * ddx1 + ddy1 * ddy1) > 1.0)
                                {
                                    // bilinear interpolation, pixval1 (radial)
                                    float xf = ii1f - ii1;
                                    float yf = jj1f - jj1;

                                    float v = data.image[IDin].array.F[kk * xysize + jj1 * xsize + ii1];
                                    float vx = data.image[IDin].array.F[kk * xysize + jj1 * xsize + ii1 + 1];
                                    float vy = data.image[IDin].array.F[kk * xysize + (jj1 + 1) * xsize + ii1];
                                    float vxy = data.image[IDin].array.F[kk * xysize + (jj1 + 1) * xsize + ii1 + 1];

                                    float pixval1 =
                                        (1.0 - xf) * (1.0 - yf) * v
                                        + (1.0 - xf) * yf * vy
                                        + xf * (1.0 - yf) * vx
                                        + xf * yf * vxy;


                                    // bilinear interpolation, pixval2 (angular)
                                    xf = ii2f - ii2;
                                    yf = jj2f - jj2;
                                    v = data.image[IDin].array.F[kk * xysize + jj2 * xsize + ii2];
                                    vx = data.image[IDin].array.F[kk * xysize + jj2 * xsize + ii2 + 1];
                                    vy = data.image[IDin].array.F[kk * xysize + (jj2 + 1) * xsize + ii2];
                                    vxy = data.image[IDin].array.F[kk * xysize + (jj2 + 1) * xsize + ii2 + 1];

                                    float pixval2 =
                                        (1.0 - xf) * (1.0 - yf) * v
                                        + (1.0 - xf) * yf * vy
                                        + xf * (1.0 - yf) * vx
                                        + xf * yf * vxy;


                                    if((pixval + pixval1 > pixfluxlim) && (pixval + pixval2 > pixfluxlim))
                                    {

                                        double pixvaldiff1 = (pixval - pixval1); // / sqrt(pixval + pixval1);
                                        double pixvaldiff2 = (pixval - pixval2); // / sqrt(pixval + pixval2);
                                        errfit1 += (pixvaldiff1 * pixvaldiff1);
                                        errfit2 += (pixvaldiff2 * pixvaldiff2);
                                        errfitcnt++;
                                    }
                                }
                            }
                        }

                        if(FASTMODE == 0)
                        {
                            data.image[IDout].array.F[9 * xysize + jjrc * xsize + iirc] = errfit1;
                            data.image[IDout].array.F[10 * xysize + jjrc * xsize + iirc] = errfit2;


                            data.image[IDout].array.F[11 * xysize + jjrc * xsize + iirc] = errfit1 -
                                    errfit2;
                            data.image[IDout].array.F[12 * xysize + jjrc * xsize + iirc] = errfitcnt;
                        }

                        float tmpval = errfit1 / errfitcnt;


                        if(tmpval < radcentval)
                        {
                            radcentx = iirc;
                            radcenty = jjrc;
                            radcentval = tmpval;
                        }


                        data.image[IDoutCent].array.F[kk * xysize + jjrc * xsize + iirc] =
                            tmpval;

                    }


                }
            }

            if(rcstep == 1)
            {
                OKloop = 0;
            }
            rcstep /= 2;

            rcradius *= 0.5;
            xcent1 = radcentx;
            ycent1 = radcenty;



            iircmin = (int)(xcent1 - rcradius);
            iircmax = (int)(xcent1 + rcradius);

            jjrcmin = (int)(ycent1 - rcradius);
            jjrcmax = (int)(ycent1 + rcradius);

            if(iircmin < 0)
            {
                iircmin = 0;
            }
            if(iircmax > xsize)
            {
                iircmax = xsize;
            }

            if(jjrcmin < 0)
            {
                jjrcmin = 0;
            }

            if(jjrcmax > ysize)
            {
                jjrcmax = ysize;
            }

        }




        free(iiarray);
        free(jjarray);
        free(iijjarray);
        free(pvarray);

        fprintf(fpout, "%5u   %8.3f %8.3f    %8.3f %8.3f   %12.3g\n", kk, photcentx,
                photcenty, radcentx, radcenty, flux);

    }


    fclose(fpout);

    return RETURN_SUCCESS;
}
