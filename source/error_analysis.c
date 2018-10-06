#include "header_rotation.h"

//probably_wrong_value_of_rotation
  double *pwvor=NULL;
  double *net_rpy;

static double getError(double *net_rpy){
  return fabs(net_rpy[ax]-pwvor[ax]);
}

static double *gpwvor(char *files[]){
  return getRingRotation(files,start_frame,end_frame,1);
}

int getBestStepSize(char *files[],int first_step_size,int last_step_size){
  int curr_step_size,best_step_size;
  double max_error=-1,error;
  FILE *out_file=fopen("output/step_size_vs_error.csv","w");

  pwvor=gpwvor(files);
  for(curr_step_size=first_step_size;curr_step_size<=last_step_size;curr_step_size++){
    //printf("At Step Size = %d\n",curr_step_size);
    net_rpy=getDirectionality(files,start_frame,end_frame,curr_step_size);
    error=getError(net_rpy);
    if(error>max_error){
      max_error=error;
      best_step_size=curr_step_size;
    }
    
    fprintf(out_file,"%d    %f    %f    %f\n",curr_step_size,pwvor[ax],net_rpy[ax],error);

    //mlfix
    free(net_rpy);
  }
  fclose(out_file);

  //mlfix
  free(pwvor);
  return best_step_size;
}
