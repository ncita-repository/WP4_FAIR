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

#define MAX_IM_COUNT NHIST-3
#define BA_PLOT_HIST NHIST-1

Imrect *all_ims[MAX_IM_COUNT];
Imrect *mean_im = NULL;
Imrect *var_im = NULL;
Imrect *roi_mask = NULL;

double *renorm_factors = NULL;
double *adj_rms = NULL;

double adj_rms_thresh = 9999999999;

int    im_count = 0;

int    skip = 1;

double hist_min = 0;
double hist_max = 100;
int    hist_bins = 50;

double ba_max = 10;
int    ba_bins = 50;

double keep_percentage = 5;

int    baseline_start = 0;
int    baseline_end = 0;

int    ox_start = 0;
int    ox_end = 0;

static void pop_mask_proc()
{
   Imrect *im;
   int type;

   if (mean_im == NULL)
   {
      format("Need to have loaded images first\n");
      return;
   }

   if (stack_check_types(IMRECT, NULL) == false)
   {
      format("No mask image on stack\n");
      return;
   }

   if (roi_mask != NULL)
   {
      im_free(roi_mask);
      roi_mask = NULL;
      format("Old mask freed\n");
   }
  
   im = (Imrect*)stack_pop(&type);
   roi_mask = im_cast(im, int_v);
   im_free(im);

   if (mean_im->region->lx != roi_mask->region->lx ||
       mean_im->region->ly != roi_mask->region->ly ||
       mean_im->region->ux != roi_mask->region->ux ||
       mean_im->region->uy != roi_mask->region->uy)
   {
      format("Rejecting mask due to region mismatch\n");
      im_free(roi_mask);
   }

   format("Mask accepted\n");
}

static void pop_images_proc()
{
   Imrect *im, *imf;
   shistogram **hists;
   int i, x, y, lx, ux, ly, uy;
   int type;
   double val;

   hists = hist_vec();

   for (i=0; i<MAX_IM_COUNT; i++)
   {
      if (all_ims[i] != NULL) im_free(all_ims[i]);
      all_ims[i] = NULL;
   }

   if (mean_im != NULL)
   {
      im_free(mean_im);
      mean_im = NULL;
   }

   i=0;
   while (i < MAX_IM_COUNT)
   {
      if (stack_check_types(IMRECT, NULL) == false)
      {
         format("Found %i images\n", i);
         format("No more images on stack\n");
         break;
      }
      
      im = (Imrect*)stack_pop(&type);
      imf = im_cast(im, float_v);
      im_free(im);
      all_ims[i] = imf;

      if (i == 0)
      {
         mean_im = im_copy(imf);
         lx = mean_im->region->lx; ly = mean_im->region->ly;
         ux = mean_im->region->ux; uy = mean_im->region->uy;
      }
      else if (lx != all_ims[i]->region->lx || ly != all_ims[i]->region->ly ||
               ux != all_ims[i]->region->ux || uy != all_ims[i]->region->uy)
      {
         format("Rejecting image due to region mismatch\n");
         im_free(all_ims[i]);
         continue;
      }

      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            val = IM_FLOAT(imf, y, x);
            if (i > 0) IM_FLOAT(mean_im, y, x) += val;
         }
      }
      i++;
   }

   im_count = i;

   for (y=ly; y<uy; y++)
   {
      for (x=lx; x<ux; x++)
      {
         val = IM_FLOAT(mean_im, y, x);
         IM_FLOAT(mean_im, y, x) = val/(double)im_count;
      }
   }

  /* stack_push(im_copy(mean_im), IMRECT, im_free);*/

  /* format("Mean image pushed to stack\n"); */
}

static void push_images_proc()
{
   int i;

   for (i=im_count-1; i>=0; i--)
   {
      if (all_ims[i] != NULL)
      {
         if (adj_rms != NULL)
         {
            if (adj_rms[i] <= adj_rms_thresh)
            {
               stack_push(im_copy(all_ims[i]), IMRECT, im_free);
               format("%i\n", i);
            }
         }
         else
         {
            stack_push(im_copy(all_ims[i]), IMRECT, im_free);
            format("%i\n", i);
         }
      }
   }

  /* stack_push(im_copy(mean_im), IMRECT, im_free);*/

   format("Images pushed\n");
}

static void hist_plot_proc()
{
   shistogram **hists;
   int i, x, y, lx, ux, ly, uy;
   int type;
   double val;

   hists = hist_vec();

   if (mean_im == NULL)
   {
      error("No images to create histograms from\n", non_fatal);
      return;
   }

   lx = mean_im->region->lx; ly = mean_im->region->ly;
   ux = mean_im->region->ux; uy = mean_im->region->uy;

   for (i=0; i<im_count; i++)
   {
      if (hists[i] != NULL) hfree(hists[i]);
      hists[i] = hbook1(i, "stack image hist", hist_min, hist_max, hist_bins);

      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               val = IM_FLOAT(all_ims[i], y, x);
               hfill1(hists[i], val, 1.0);
            }
         }
      }

      format("%i mean %f\n", i, hists[i]->mean);

   }

   format("Histograms populated\n");
}

static void ba_from_mean_proc()
{
   int lx, ly, ux, uy, x, y;
   double val, ref;
   int i;
   shistogram **hists;

   if (mean_im == NULL)
   {
      error("No mean image to create plot from\n", non_fatal);
      return;
   }

   hists = hist_vec();

   lx = mean_im->region->lx; ly = mean_im->region->ly;
   ux = mean_im->region->ux; uy = mean_im->region->uy;

   if (hists[BA_PLOT_HIST] != NULL) hfree(hists[BA_PLOT_HIST]);
   hists[BA_PLOT_HIST] = hbook2(BA_PLOT_HIST, "BA plot", hist_min, hist_max, hist_bins, -ba_max, ba_max, ba_bins);

   for (i=0; i<im_count; i++)
   {
      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               val = IM_FLOAT(all_ims[i], y, x);
               ref = IM_FLOAT(mean_im, y, x);
               hfill2(hists[BA_PLOT_HIST], ref, val-ref, 1.0);
            }
         }
      }
   }

   format("Done\n");
}

static void ba_from_interp_proc()
{
   format("THIS DOES NOT DO ANYTHING RIGHT NOW!\n");
}

static void var_from_mean_proc()
{
   int lx, ly, ux, uy, x, y;
   double val, ref;
   int i;

   if (mean_im == NULL)
   {
      error("No mean image to create plot from\n", non_fatal);
      return;
   }

   if (var_im != NULL)
   {
      im_free(var_im);
      var_im = NULL;
   }

   var_im = im_alloc(mean_im->height, mean_im->width, mean_im->region, mean_im->vtype);

   lx = mean_im->region->lx; ly = mean_im->region->ly;
   ux = mean_im->region->ux; uy = mean_im->region->uy;

   for (i=0; i<im_count; i++)
   {
      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               val = IM_FLOAT(all_ims[i], y, x);
               ref = IM_FLOAT(mean_im, y, x);
               IM_FLOAT(var_im, y, x) += ((val-ref)*(val-ref))/(double)im_count;
            }
         }
      }
   }

   stack_push(im_copy(var_im), IMRECT, im_free);
   format("Variance image pushed to stack\n");
}

static void var_from_interp_proc()
{
   int lx, ly, ux, uy, x, y;
   double val, ref;
   int i;

   if (mean_im == NULL)
   {
      error("No mean image to create plot from\n", non_fatal);
      return;
   }

   if (var_im != NULL)
   {
      im_free(var_im);
      var_im = NULL;
   }

   var_im = im_alloc(mean_im->height, mean_im->width, mean_im->region, mean_im->vtype);

   lx = mean_im->region->lx; ly = mean_im->region->ly;
   ux = mean_im->region->ux; uy = mean_im->region->uy;

   if (im_count >= 3) for (i=1; i<im_count-1; i++)
   {
      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               ref = IM_FLOAT(all_ims[i-1], y, x);
               val = IM_FLOAT(all_ims[i], y, x);
               ref += IM_FLOAT(all_ims[i+1], y, x);
               ref /= 2.0;
               IM_FLOAT(var_im, y, x) += ((val-ref)*(val-ref))/(double)(im_count-1);
            }
         }
      }
   }

   stack_push(im_copy(var_im), IMRECT, im_free);
   format("Variance image pushed to stack\n");
}

static Imrect *mean_image(int start, int end)
{
   Imrect *im;
   int i, lx, ux, ly, uy, x, y, n;
   double val;

   n = (end-start) + 1;

   for (i=start; i<=end; i++)
   {
      if (i == start)
      {
         im = im_times(1.0/(double)n, all_ims[i]);
         lx = im->region->lx; ly = im->region->ly;
         ux = im->region->ux; uy = im->region->uy;
      }
      else for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            val = IM_FLOAT(all_ims[i], y, x);
            IM_FLOAT(im, y, x) += val/(double)n;
         }
      }
   }

   if (roi_mask != NULL) for (y=ly; y<uy; y++)
   {
      for (x=lx; x<ux; x++)
      {
         if (IM_INT(roi_mask, y, x) == 0) IM_FLOAT(im, y, x)  = 0;
      }
   }

   return im;
}

static void mean_baseline_proc()
{
   stack_push(mean_image(baseline_start, baseline_end), IMRECT, im_free);
}

static void mean_ox_proc()
{
   stack_push(mean_image(ox_start, ox_end), IMRECT, im_free);
}

static double percentile_location(Imrect *ref, Imrect *im)
{
   Imrect *diff;
   shistogram *diff_hist;
   double val, total, thresh;
   int lx, ux, ly, uy, x, y;

   diff = im_diff(ref, im);

   val = im_locate_max(diff, &x, &y); /* could find outliers, not good! */
   diff_hist = hbook1(0, "diff", 0, val/2.0, 100);

   lx = diff->region->lx; ly = diff->region->ly;
   ux = diff->region->ux; uy = diff->region->uy;

   for (y=ly; y<uy; y++)
   {
      for (x=lx; x<ux; x++)
      {
         if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
         {
            val = IM_FLOAT(diff, y, x);
            hfill1(diff_hist, abs(val), 1.0);
         }
      }
   }

   hintegf(diff_hist);

   total = diff_hist->array[0][diff_hist->xbins-1];
   for (x=0; x<diff_hist->xbins; x++)
   {
      val = diff_hist->array[0][x];
      if (val/total*100 >= keep_percentage) break;
   }

   thresh = diff_hist->xincr*x;
   format("kp at %i, %f, hist max %f\n", x, thresh, diff_hist->xmax);

   im_free(diff);

/*   im_free(baseline);*/



hist_vec()[0] = diff_hist;
/* hfree(diff_hist); */
}








static void check_norm_proc()
{
   double profile, smoothed_profile, val;
   int lx, ux, ly, uy, x, y, i, j, n;

   if (renorm_factors != NULL)
   {
      dvector_free(renorm_factors, 0);
      renorm_factors = NULL;
   }
   renorm_factors = dvector_alloc(0, im_count);

   if (adj_rms != NULL)
   {
      dvector_free(adj_rms, 0);
      adj_rms = NULL;
   }
   adj_rms = dvector_alloc(0, im_count);

   lx = all_ims[0]->region->lx; ly = all_ims[0]->region->ly;
   ux = all_ims[0]->region->ux; uy = all_ims[0]->region->uy;

   n = 0;
   for (y=ly; y<uy; y++) for (x=lx; x<ux; x++)
   {
      if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0) n++;
   }

   format("t, mean profile, smoothed profile, renorm factor, adj rms\n");
   printf("t, mean profile, smoothed profile, renorm factor, adj rms\n");
   
   for (j=0; j<im_count; j++)
   {
      profile = 0.0; smoothed_profile = 0.0;
      i = j;
      if (j < 2) i = 2;
      if (j > im_count-3) i = im_count-3;

      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               val = IM_FLOAT(all_ims[i-2], y, x);
               smoothed_profile += val/(double)(n*4);

               val = IM_FLOAT(all_ims[i-1], y, x);
               smoothed_profile += val/(double)(n*4);

               val = IM_FLOAT(all_ims[i], y, x);
               smoothed_profile += val/(double)(n*4);

               val = IM_FLOAT(all_ims[i+1], y, x);
               smoothed_profile += val/(double)(n*4);

               val = IM_FLOAT(all_ims[i+2], y, x);
               smoothed_profile += val/(double)(n*4);
            
               val = IM_FLOAT(all_ims[j], y, x);    
               profile += val/(double)n;
               smoothed_profile -= val/(double)(n*4);

               if (j>0 && j<im_count-1)
               {
                  val = (IM_FLOAT(all_ims[j-1], y, x)+IM_FLOAT(all_ims[j+1], y, x))/2;
                  val -= IM_FLOAT(all_ims[j], y, x);
                  adj_rms[j] += val*val;
               }
            }
         }
      }

      adj_rms[j] = sqrt(adj_rms[j]/n);
      renorm_factors[j] = smoothed_profile-profile;

      format("%i, %f, %f, %f, %f\n", j, profile, smoothed_profile, renorm_factors[j], adj_rms[j]);
      printf("%i, %f, %f, %f, %f\n", j, profile, smoothed_profile, renorm_factors[j], adj_rms[j]);
   }
}

static void renorm_proc()
{
   int lx, ux, ly, uy, x, y, i;
   double val;

   if (renorm_factors == NULL) return;

   lx = all_ims[0]->region->lx; ly = all_ims[0]->region->ly;
   ux = all_ims[0]->region->ux; uy = all_ims[0]->region->uy;

   for (i=0; i<im_count; i++)
   {
      for (y=ly; y<uy; y++)
      {
         for (x=lx; x<ux; x++)
         {
            if (roi_mask == NULL || IM_INT(roi_mask, y, x) != 0)
            {
               val = IM_FLOAT(all_ims[i], y, x);
               IM_FLOAT(all_ims[i], y, x) = val+renorm_factors[i]; 
            }
         }
      }
   }

   format("Renormed\n");
}

static void z_score_proc()
{
   Imrect *baseline, *ox;
   Imrect *bl_var, *ox_var;
   Imrect *diff, *diff_var, *diff_sd;
   Imrect *zscores;

   if (var_im == NULL)
   {
      error("Need to compute var from interp image first\n", non_fatal);
      return;
   }

   baseline = mean_image(baseline_start, baseline_end);
   ox = mean_image(ox_start, ox_end);
   diff = im_diff(ox, baseline);

   bl_var = im_times(1.0/(1.5*(1+baseline_end-baseline_start)), var_im);
   ox_var = im_times(1.0/(1.5*(1+ox_end-ox_start)), var_im);
   diff_var = im_sum(bl_var, ox_var);
   diff_sd = im_sqrt(diff_var);

/*   zscores = */

   stack_push(diff_sd, IMRECT, im_free);
   stack_push(diff, IMRECT, im_free);
/*   stack_push(zscores, IMRECT, im_free);*/

   format("Pushed diff sd, diff and zscores\n");

   im_free(diff_var);
   im_free(bl_var);
   im_free(ox_var);
   im_free(baseline);
   im_free(ox);
}

static void save_stim_proc()
{
  Sequence *seq = seq_get_current();
  float    *time, *flow;
  int       imptrlx, imptrux, imptrly, imptruy, imptrlz, imptruz;
  int       length,  k;
  float    *stim;

  if ((stim=get_stimulus()) == NULL)
    return;

  seq_slice_init(seq);
  seq_limits(&imptrlz, &imptruz, &imptrly, &imptruy, &imptrlx, &imptrux);
  
  length = imptruz-imptrlz;    
  flow = (float *)fvector_alloc(0, length);
  time = (float *)fvector_alloc(0, length); 

  printf("yes, it's doing something\n");

  for (k = 0; k < length; k++)
    {
      if (!isnan((double)stim[k]))
        flow[k] = stim[k];
      else
        flow[k] = 0.0;
      time[k] = (float)k;
      fprintf(stdout," %f %f \n", time[k], flow[k]);
    }

   fvector_free(flow, 0);
   fvector_free(time, 0);
}

/********** Tool creation **********/

void            image_info_tool(int x, int y)
{
    int i;
    static void *tool = NULL;

    if (tool)
    {
	tw_show_tool(tool);
	return;
    }

    for (i=0; i<MAX_IM_COUNT; i++) all_ims[i] = NULL;

    tool = (void *)tw_tool("Image Info Tool", x, y);

    tw_fglobal("min:", &hist_min, 4);
    tw_fglobal("max:", &hist_max, 4);
    tw_iglobal("bins:", &hist_bins, 4);
    tw_iglobal("skip:", &skip, 4);
    tw_newrow();

    tw_fglobal("ba max:", &ba_max, 4);
    tw_iglobal("ba bins:", &ba_bins, 4);
    tw_fglobal("keep %:", &keep_percentage, 4);
    tw_newrow();

    tw_iglobal("baseline start:", &baseline_start, 4);
    tw_iglobal("baseline end:", &baseline_end, 4);
    tw_newrow();

    tw_iglobal("ox start:", &ox_start, 4);
    tw_iglobal("ox end:", &ox_end, 4);
    tw_newrow();

    tw_fglobal("Adj rms thresh:", &adj_rms_thresh, 8);
    tw_newrow();
        
    tw_button("Pop images", pop_images_proc, NULL);
    tw_button("Push images", push_images_proc, NULL);
    tw_button("Pop mask", pop_mask_proc, NULL);
    tw_button("check norm", check_norm_proc, NULL);
    tw_button("Renorm", renorm_proc, NULL);
    tw_newrow();

    tw_button("Histograms", hist_plot_proc, NULL);
    tw_button("BA from mean image", ba_from_mean_proc, NULL);
    tw_button("Var from mean image", var_from_mean_proc, NULL);
    tw_newrow();

    tw_button("BA from interp", ba_from_interp_proc, NULL);
    tw_button("Var from interp", var_from_interp_proc, NULL);
    tw_button("Mean baseline", mean_baseline_proc, NULL);
    tw_button("Mean ox", mean_ox_proc, NULL);

    tw_button("Z-score", z_score_proc, NULL);

    tw_button("Save stim", save_stim_proc, NULL);


    tw_end_tool();
}
