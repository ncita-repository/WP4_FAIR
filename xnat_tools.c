#include <stdio.h>
#include <sys/param.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <tina/sys/sysDef.h>
#include <tina/sys/sysPro.h>
#include <tina/math/mathDef.h>
#include <tina/math/mathPro.h>
#include <tina/image/imgDef.h>
#include <tina/image/imgPro.h>
#include <tina/geometry/geomDef.h>
#include <tina/geometry/geomPro.h>
#include <tinatool/draw/drawDef.h>
#include <tinatool/draw/drawPro.h>
#include <tinatool/draw/drawDef.h>
#include <tinatool/draw/drawPro.h>
#include <tinatool/wdgts/wdgtsDef.h>
#include <tinatool/wdgts/wdgtsPro.h>
#include <tinatool/tlbase/tlbase_InfrDef.h>
#include <tinatool/tlbase/tlbase_InfrPro.h>
#include <tinatool/tlvision/tlvisEdge_epi_mouse.h>

#define NOISE 0
#define MEAN  1
#define VAR   2
#define RANGE 3

/*********** private ***********/

static char filename[1024];

/************ funcs ************/

static void record_output(int field, float val)
{
/* LONG-TERM PLAN TO SAVE IN JSON FIELD FORMAT??? */

   FILE *fp;

   if ((fp = fopen(filename, "at")) == NULL)
   {
      error("Unable to open file\n", non_fatal);
      return;
   }

   switch (field)
   {
      case NOISE: fprintf(fp, "NOISE, "); break;
      case MEAN:  fprintf(fp, "MEAN, "); break;
      case VAR:   fprintf(fp, "VAR, "); break;
      case RANGE: fprintf(fp, "RANGE, "); break;
   }

   fprintf(fp, "%f\n", val);

   fclose(fp);
}

static void get_image_noise()
{
   Imrect *im;
   double sigma_x, sigma_y, sigma;
   int type;

   if (stack_check_types(IMRECT, NULL) == false)
   {
      error("Wrong type on stack", warning);
      return;
   }

   im = (Imrect *) stack_pop(&type);

   sigma_x = imf_diffx_noise(im, im->region);
   sigma_y = imf_diffy_noise(im, im->region);
   sigma = (sigma_x + sigma_y)/2;

   stack_push(im, IMRECT, im_free);

   format("sigma x: %f, sigma y: %f, sigma: %f\n", sigma_x, sigma_y, sigma);

   record_output(NOISE, sigma);
}

static void get_image_var()
{
   int type;
   int lx, ly, ux, uy, x, y;
   Imrect *im, *imf;

   double mean, m2, val, delta1, delta2, var;
   int n;

   if (stack_check_types(IMRECT, NULL) == false)
   {
      format("No image on stack\n");
      return;
   }
      
   im = (Imrect*)stack_pop(&type);
   imf = im_cast(im, float_v);

   lx = im->region->lx; ly = im->region->ly;
   ux = im->region->ux; uy = im->region->uy;

   n = 0; mean = 0.0; m2 = 0.0;

   for (y=ly; y<uy; y++)
   {
      for (x=lx; x<ux; x++)
      {
         val = IM_FLOAT(imf, y, x);
         n++;
         delta1 = val - mean;
         mean += delta1 / n;
         delta2 = val - mean;
         m2 += delta1 * delta2;
      }
   }

   var = m2 / (n-1);
   
   im_free(imf);
   stack_push(im, IMRECT, im_free);

   format("mean: %f, var: %f\n", mean, var);

   record_output(MEAN, mean);
   record_output(VAR, var);
}

static void get_image_range()
{
   int type;
   Imrect *im, *imf;
   float min, max;
   double range;

   if (stack_check_types(IMRECT, NULL) == false)
   {
      format("No image on stack\n");
      return;
   }
      
   im = (Imrect*)stack_pop(&type);
   imf = im_cast(im, float_v);

   imf_minmax(imf, &min, &max);
   range = max-min;

   im_free(imf);
   stack_push(im, IMRECT, im_free);

   format("min: %f, max %f, range %f\n", min, max, range);

   record_output(RANGE, range);
}

/********** Tool creation **********/

void xnat_tools(int x, int y)
{
    static void *tool = NULL;

    if (tool)
    {
	tw_show_tool(tool);
	return;
    }

    tool = (void *)tw_tool("XNAT Tools", x, y);

    tw_sglobal("Output file:", filename, 42);
    tw_newrow();

    tw_button("noise", get_image_noise, NULL);
    tw_button("variance", get_image_var, NULL);
    tw_button("range", get_image_range, NULL);

    tw_end_tool();
}
