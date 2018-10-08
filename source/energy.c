#include "header_rotation.h"


double getFrameRangeEnergy(char *files[],int start_frame,int end_frame,int step_size){
  int i=0,j=0,f1=-1,f2=-1,is_read=-1;
  FILE *file1,*file2;
  //FILE *mix_out=fopen("output/mix_out.dat","w");
  //FILE *pos_out=fopen("output/pos_out.dat","w");
  //FILE *neg_out=fopen("output/neg_out.dat","w");
  //FILE *torque_out=fopen("output/torque.dat","w");
  double rot_ke,trans_ke,**tmp_cords;//*net_rpy,*rpy,*net_torque,*trpy;

  //net_rpy=(double *)malloc(sizeof(double)*3);
  //for(i=0;i<3;i++)net_rpy[i]=0;
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

  for(f2=f1+step_size;f2<=end_frame;f2+=step_size){
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

    //printf("First\n");
    //printCords(f2cords,c);

    printf("\nstep_size = %d frame = %d\n",step_size,f1);
    rot_ke=getRingRotKE(files,f1,f2,0);
    trans_ke=getRingTransKE(files,f1,f2,0);
    printf("rot_ke=%e trans_ke=%e\n",rot_ke,trans_ke);
    /*rpy=getRingRotation(files,f1,f2,0);
    for(i=0;i<3;i++)printf("%f ",rpy[i]);printf("\n");

    //printf("Second\n");
    //printCords(f2cords,c);

    //track_rotation
    if(0){
    trpy=getTrackRotation(files,f1,f2,track_start_atom_no,track_end_atom_no,0);
    for(i=0;i<3;i++)printf("%f ",trpy[i]);printf("\n");
    for(i=0;i<3;i++)rpy[i]-=trpy[i];
    for(i=0;i<3;i++)printf("%f ",rpy[i]);printf("\n");

    //ml_fix
    free(trpy);
    }

    //printf("Third\n");
    //printCords(f2cords,c);

    //torque
   
    net_torque=getTorqueInFrame(f2);
    for(i=0;i<3;i++)printf("%f ",net_torque[i]);printf("\n");
    writeTorqueOutput(torque_out,f2,net_torque,rpy[ax]);
    

    writeOutput(mix_out,f1,rpy);
    if(rpy[ax]<0)writeOutput(neg_out,f1,rpy);
    else writeOutput(pos_out,f1,rpy);
    for(i=0;i<3;i++) net_rpy[i]+=rpy[i];*/
    f1=f2;
    file1=file2;
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        f1cords[i][j]=tmp_cords[i][j];
      }
    }

    //mlfix
    //free(rpy);
    //free(trpy);
  }
  fclose(file2);
  //fclose(mix_out);
  //fclose(pos_out);
  //fclose(neg_out);

  //mlfix
  for(i=0;i<atoms;i++)free(tmp_cords[i]);
  free(tmp_cords);
  //free(rpy);

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
  //printf("MI=%e w=%e\n",getMI(n),w);
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

/*
double getTrans(double v1[],double v2[]){
  int i=0;
  double d=0.0;
  double trans;*trans=(double *)malloc(sizeof(double)*3);
  trans=getMag(v2)-getMag(v1);
  trans=getuv1(v2);
  for(i=0;i<3;i++){
    trans[i]*=d;
  }
  return trans;
}
*/

double getAtmTransKE(int n){
  int i=0;
  double trans,v;
  printf("f1 ");for(i=0;i<3;i++)printf("%f ",f1cords[n][i]);printf("\n");
  printf("f2 ");for(i=0;i<3;i++)printf("%f ",f2cords[n][i]);printf("\n");
  trans=getTrans(f1cords[n],f2cords[n]);
  v=trans*angstorm/(timeStep*femto);
  printf("t=%e\n",trans);
  return 0.5*getAtmass(n,c)*amu*pow(v,2);

}


double getRingTransKE(char **files,int f1,int f2,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,ke=0;
  /*if(test){
    printCords(f1cords,c);
    printCords(f2cords,c);
  }*/
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
  if(!test){
    printCords(f1cords,c);
    printCords(f2cords,c);
  }
  for(i=0;i<ring_atom;i++){
    ke=+getAtmTransKE(i);
  }
  
  //mlfix
  free(f1com);
  free(f2com);

  return ke;
}

