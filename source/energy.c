#include "header_rotation.h"

/*void getEnergyStats(double *rot_ke,double *trans_ke){
  int i=0;
  double d=0,avg_rot_ke=0,avg_trans_ke=0,effi[ring_atom],avg_effi=0;
  for(i=0;i<ring_atom;i++){
    avg_rot_ke+=rot_ke[i];
    avg_trans_ke+=trans_ke[i];
    effi[i]=rot_ke[i]*100/(rot_ke[i]+trans_ke[i]);
    avg_effi+=effi[i];
  }
  
  printf("\navg rot=%e trans=%e effi=%f\n",avg_rot_ke/ring_atom,avg_trans_ke/ring_atom,avg_effi/ring_atom);
}*/

double **createData(double *a,double *b,double *c,double *d,int n){
  int i=0;
  double **data;
  data=join1dArrays_d(a,b,n);
  data=join21dArrays_d(data,c,n,2);
  data=join21dArrays_d(data,d,n,3);
  return data;
}

double getFrameRangeEnergy(char *files[],int start_frame,int end_frame,int step_size){
  int i=0,j=0,f1=-1,f2=-1,is_read=-1,itr=1,data_points=(end_frame-start_frame)/step_size;
  FILE *file1,*file2;
  //FILE *mix_out=fopen("output/mix_out.dat","w");
  double **data,frame_buff[data_points],rot_ke_buff[data_points],trans_ke_buff[data_points],effi_buff[data_points],rot_ke=0,trans_ke=0,effi=0,avg_rot_ke=0,avg_trans_ke=0,avg_effi=0,**tmp_cords;

  tmp_cords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
  *(tmp_cords+i)=(double *)malloc(sizeof(double)*3);
  }
  f1=start_frame;
  file2=file1=fopen(f1<=file_limit?files[0]:files[1],"r");
  is_read=readFile(file1,f1cords,c,f1);
  if(!is_read){
    printf("Error in reading files at first frame\n");
    return -1;
  }else{
    //printf("Found frame %d\n",f1);
  }

  for(f2=f1+step_size,itr=1;f2<=end_frame;f2+=step_size,itr++){
    if((f1<=file_limit?files[0]:files[1])!=(f2<=file_limit?files[0]:files[1])){
      file2=fopen(f2<=file_limit?files[0]:files[1],"r");
      //printf("file2 pointer changed\n");
      fclose(file1);
    }else{
      file2=file1;
      //printf("file2 pointer is same\n");
    }
    is_read=readFile(file2,f2cords,c,f2);
    if(!is_read){
      printf("Error in reading files....reached till frame %d \n",f1);
      return -1;
    }else{
      //printf("Found frame %d\n",f2);
    }
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        tmp_cords[i][j]=f2cords[i][j];
      }
    }

    printf("\nstep_size = %d frame = %d\n",step_size,f1);
    rot_ke=getRingRotKE(files,f1,f2,0);
    trans_ke=getRingTransKE(files,f1,f2,0);
    effi=rot_ke*100/(rot_ke+trans_ke);
    printf("rot_ke=%e trans_ke=%e effi=%f\n",rot_ke,trans_ke,effi);
    avg_rot_ke+=rot_ke;
    avg_trans_ke+=trans_ke;
    avg_effi+=effi;
   
    //buffers
    rot_ke_buff[itr-1]=rot_ke;
    trans_ke_buff[itr-1]=trans_ke;
    effi_buff[itr-1]=effi;
    frame_buff[itr-1]=f1;

    f1=f2;
    file1=file2;
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        f1cords[i][j]=tmp_cords[i][j];
      }
    }
   
  }
  itr--;
  printf("\n%d avg rot=%e trans=%e effi=%f\n",itr,avg_rot_ke/itr,avg_trans_ke/itr,avg_effi/itr);
  data=createData(frame_buff,rot_ke_buff,trans_ke_buff,effi_buff,data_points);
  writeCsv("output/test_csv.csv","frame,rot_ke,trans_ke,effi",data,data_points,4);
  //getEnergyStats(rot_ke,trans_ke);
  fclose(file2);

  //mlfix
  for(i=0;i<atoms;i++)free(tmp_cords[i]);
  free(tmp_cords);

  return 1;
}


double getMI(int n){
  int i;
  double d,I;
  return getAtmass(n,c)*amu*pow(getMag(f1cords[n])*angstorm,2);  
}


double getAtmRotKE(int n){
  double theta,w;

  theta=getAngleR(f1cords[n],f2cords[n]);
  w=theta/(timeStep*femto);
  return 0.5*getMI(n)*pow(w,2);
  
}

double getRingRotKE(char **files,int f1,int f2,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,ke=0;
  if(test){
    FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
    FILE *file2=fopen(f2<=file_limit?files[0]:files[1],"r");
    isf1=readFile(file1,f1cords,c,f1);
    isf2=readFile(file2,f2cords,c,f2);
    if(!isf1||!isf2){
      printf("Error in reading files\n");
      return -1;
    }else{
      printf("Found frame %d and %d\n",f1,f2);
    }
    fclose(file1);
    fclose(file2);

    //debug
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        orig_f1cords[i][j]=f1cords[i][j];
      }
    }
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        orig_f2cords[i][j]=f2cords[i][j];
      }
    }
  }

  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++)f2com[i]=f2com[i]-f1com[i];
  shiftOrigin(f1com,f2com);

  for(i=0;i<ring_atom;i++){
    ke=+getAtmRotKE(i);
  }

  //mlfix
  free(f1com);
  free(f2com);
  
  return ke;
}


double getAtmTransKE(int n){
  int i=0;
  double trans,v;
  trans=getTrans(f1cords[n],f2cords[n]);
  v=trans*angstorm/(timeStep*femto);
  return 0.5*getAtmass(n,c)*amu*pow(v,2);
}


double getRingTransKE(char **files,int f1,int f2,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,ke=0;
  
  if(test){
    FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
    FILE *file2=fopen(f2<=file_limit?files[0]:files[1],"r");
    isf1=readFile(file1,f1cords,c,f1);
    isf2=readFile(file2,f2cords,c,f2);
    if(!isf1||!isf2){
      printf("Error in reading files\n");
      return -1;
    }else{
      printf("Found frame %d and %d\n",f1,f2);
    }
    fclose(file1);
    fclose(file2);

    //debug
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        orig_f1cords[i][j]=f1cords[i][j];
      }
    }
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        orig_f2cords[i][j]=f2cords[i][j];
      }
    }
  
    f1com=getCom(f1cords,ring_atom);
    f2com=getCom(f2cords,ring_atom);
    for(i=0;i<3;i++)f2com[i]=f2com[i]-f1com[i];
    shiftOrigin(f1com,f2com);
  }

  for(i=0;i<ring_atom;i++){
    ke=+getAtmTransKE(i);
  }
  
  //mlfix
  free(f1com);
  free(f2com);

  return ke;
}

