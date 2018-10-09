#include "header_rotation.h"

double getNetTranslation(char *files[],int start_frame,int end_frame,int step_size){  
  int i=0,j=0,f1=-1,f2=-1,is_read=-1;
  FILE *file1,*file2;
  double net_disp,**tmp_cords,disp;
  //net_disp=(double *)malloc(sizeof(double)*4);
  //for(i=0;i<3;i++)net_disp[i]=0;
  tmp_cords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
  *(tmp_cords+i)=(double *)malloc(sizeof(double)*3);
  }
  f1=start_frame;
  file1=fopen(f1<=file_limit?files[0]:files[1],"r");
  is_read=readFile(file1,f1cords,c,f1);
  if(!is_read){
    printf("Error in reading files\n");
    return -1;
  }else{
    //printf("Found frame %d\n",f1);
  }
  for(f2=f1+step_size;f2<=end_frame;f2+=step_size){
    if((f1<=file_limit?files[0]:files[1])!=(f2<=file_limit?files[0]:files[1])){
      file2=fopen(f2<=file_limit?files[0]:files[1],"r");
      //printf("file2 pointer changed\n");
      fclose(file1);
    }else{
      file2=file1;
      //printf("file2 pointer is same\n");
    }
    //file2=fopen(f2<=file_limit?files[0]:files[1],"r");
    is_read=readFile(file2,f2cords,c,f2);
    if(!is_read){
      printf("Error in reading files\n");
      return -1;
    }else{
      //printf("Found frame %d\n",f2);
    }

    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        tmp_cords[i][j]=f2cords[i][j];
      }
    }
    disp=getRingTranslation(files,f1,f2,0);
    net_disp+=disp;
    f1=f2;
    file1=file2;
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        f1cords[i][j]=tmp_cords[i][j];
      }
    }
  }
  fclose(file2);
  return net_disp;
}

double getRingTranslation(char **files,int f1,int f2,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,v,avg_disp,disp;
  //disp=(double *)malloc(sizeof(double)*3);
  //for(i=0;i<3;i++)disp[i]=0;
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
  }

  /*
  printf("\nf1 Cords %d\n",f1);
  printCords(f1cords,c);
  printf("\nf2 Cords %d\n",f2);
  printCords(f2cords,c);
  */
 
  /*
  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++)f2com[i]=f2com[i]-f1com[i];
  shiftOrigin(f1com,f2com);
  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  //createFilez(ring_atom,1,f1cords,f1com,f2com);
  //createFilez(ring_atom,2,f2cords,f1com,f2com);
  */
  //for(i=0;i<3;i++)disp[i]=f2com[i]-f1com[i];
  
  for(i=0;i<ring_atom;i++){
    v=getTrans(f1cords[i],f2cords[i]);
    printf("%f ",v);  //*
    avg_disp+=v/ring_atom;
  }printf("\n"); //*

  
    
  printf("%f ",avg_disp);printf("\n");
  return avg_disp;
}

double *getTrasformationTrans(double v1[],double v2[]){
  int i=0;
  double d=0.0;
  double *trans=(double *)malloc(sizeof(double)*3);
  d=getMag(v2)-getMag(v1);
  trans=getuv1(v2);
  for(i=0;i<3;i++){
    trans[i]*=d;
  }
  return trans;
}


//**Distance only**//

double getTrans(double v1[],double v2[]){
  int i=0;
  double *trans=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++){
    trans[i]=v2[i]-v1[i];
  }
  return getMag(trans);
}

//trans with dm == trans with COM

//double *getCOMTrans(double *com1,double *

