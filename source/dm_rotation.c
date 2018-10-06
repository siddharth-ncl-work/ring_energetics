#include "header_rotation.h"

double *dmGetDirectionality(char *files[],int start_frame,int end_frame,int step_size){
  int i=0,j=0,f1=-1,f2=-1,is_read=-1;
  FILE *file1,*file2;
  FILE *mix_out=fopen("output/mix_out.dat","w");
  FILE *pos_out=fopen("output/pos_out.dat","w");
  FILE *neg_out=fopen("output/neg_out.dat","w");
  FILE *torque_out=fopen("output/torque.dat","w");
  double *net_rpy,**tmp_cords,*rpy,*net_torque,*trpy;

  net_rpy=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++)net_rpy[i]=0;
  tmp_cords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
  *(tmp_cords+i)=(double *)malloc(sizeof(double)*3);
  }
  f1=start_frame;
  file2=file1=fopen(f1<=file_limit?files[0]:files[1],"r");
  is_read=readFile(file1,f1cords,c,f1);
  if(!is_read){
    printf("Error in reading files at first frame\n");
    return NULL;
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
      return NULL;
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
    rpy=dmGetRingRotation(files,f1,f2,0);
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
    /*
    net_torque=getTorqueInFrame(f2);
    for(i=0;i<3;i++)printf("%f ",net_torque[i]);printf("\n");
    writeTorqueOutput(torque_out,f2,net_torque,rpy[ax]);
    */

    writeOutput(mix_out,f1,rpy);
    if(rpy[ax]<0)writeOutput(neg_out,f1,rpy);
    else writeOutput(pos_out,f1,rpy);
    for(i=0;i<3;i++) net_rpy[i]+=rpy[i];
    f1=f2;
    file1=file2;
    for(i=0;i<atoms;i++){
      for(j=0;j<3;j++){
        f1cords[i][j]=tmp_cords[i][j];
      }
    }

    //mlfix
    free(rpy);
    //free(trpy);
  }
  fclose(file2);
  fclose(mix_out);
  fclose(pos_out);
  fclose(neg_out);

  //mkfix
  for(i=0;i<atoms;i++)free(tmp_cords[i]);
  free(tmp_cords);
  //free(rpy);

  return net_rpy;
}

double *dmGetRingRotation(char **files,int f1,int f2,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,*avg_rpy,*v,*dm1,*dm2;//,rpy[3]={0,0,0};
  //avg_rpy=(double *)malloc(sizeof(double)*3);
  //for(i=0;i<3;i++)avg_rpy[i]=0;
  if(test){
    FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
    FILE *file2=fopen(f2<=file_limit?files[0]:files[1],"r");
    isf1=readFile(file1,f1cords,c,f1);
    isf2=readFile(file2,f2cords,c,f2);
    if(!isf1||!isf2){
      printf("Error in reading files\n");
      return NULL;
    }else{
      printf("Found frame %d and %d\n",f1,f2);
    }
    fclose(file1);
    fclose(file2);
  }
 
  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++)f2com[i]=f2com[i]-f1com[i];
  shiftOrigin(f1com,f2com);

  //mlfix
  free(f1com);
  free(f2com);

  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  //createFilez(ring_atom,1,f1cords,f1com,f2com);
  //createFilez(ring_atom,2,f2cords,f1com,f2com);
  
  
  dm1=getDm(f1cords,charge,0,ring_atom-1,f1com);
  dm2=getDm(f2cords,charge,0,ring_atom-1,f2com);

  avg_rpy=getRPY(dm1,dm2);

  //mlfix
  free(f1com);
  free(f2com);
  free(dm1);
  free(dm2);
  for(i=0;i<3;i++)printf("%f ",avg_rpy[i]);printf("\n");
  return avg_rpy;
}




//track rotation

double *dmGetTrackRotation(char **files,int f1,int f2,int start_atom_no,int end_atom_no,int test){
  int i=0,j=0,isf1=0,isf2=0;
  double *f1com,*f2com,*avg_rpy=NULL,*v,*dm1,*dm2;//,rpy[3]={0,0,0};
  //avg_rpy=(double *)malloc(sizeof(double)*3);
  //for(i=0;i<3;i++)avg_rpy[i]=0;
  if(test){
    FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
    FILE *file2=fopen(f2<=file_limit?files[0]:files[1],"r");
    isf1=readFile(file1,f1cords,c,f1);
    isf2=readFile(file2,f2cords,c,f2);
    if(!isf1||!isf2){
      printf("Error in reading files\n");
      return NULL;
    }else{
      printf("Found frame %d and %d\n",f1,f2);
    }
    fclose(file1);
    fclose(file2);

    f1com=getCom(f1cords,ring_atom);
    f2com=getCom(f2cords,ring_atom);
    for(i=0;i<3;i++)f2com[i]=f2com[i]-f1com[i];
    shiftOrigin(f1com,f2com);

    //mlfix
    free(f1com);
    free(f2com);

    f1com=getCom(f1cords,ring_atom);
    f2com=getCom(f2cords,ring_atom);
    //createFilez(ring_atom,1,f1cords,f1com,f2com);
    //createFilez(ring_atom,2,f2cords,f1com,f2com);
  }

  dm1=getDm(f1cords,charge,track_start_atom_no,track_end_atom_no,f1com);
  dm2=getDm(f2cords,charge,track_start_atom_no,track_end_atom_no,f2com);

  avg_rpy=getRPY(dm1,dm2);
  

  //mlfix
  free(dm1);
  free(dm2);
  //free(f1com);
  //free(f2com);
  //free(v);
  //for(i=0;i<3;i++)printf("%f ",avg_rpy[i]);printf("\n");
  return avg_rpy;
}

  
/*


double **getR(double *s,double t){
  int i,j;
  double **R;
  double vv;
  R=(double **)malloc(sizeof(double*)*3);
  for(i=0;i<3;i++){
    *(R+i)=(double *)malloc(sizeof(double)*3);
  }

  vv=(1-cos(t));
  R[0][0]=s[0]*s[0]*vv+cos(t);
  R[0][1]=s[0]*s[1]*vv-s[2]*sin(t);
  R[0][2]=s[0]*s[2]*vv+s[1]*sin(t);
  R[1][0]=s[0]*s[1]*vv+s[2]*sin(t);
  R[1][1]=s[1]*s[1]*vv+cos(t);
  R[1][2]=s[1]*s[2]*vv-s[0]*sin(t);
  R[2][0]=s[0]*s[2]*vv-s[1]*sin(t);
  R[2][1]=s[1]*s[2]*vv+s[0]*sin(t);
  R[2][2]=s[2]*s[2]*vv+cos(t);

  return R;
}

double *getRPY(double *v1,double *v2){
  int i=0,j=0;
  double *s,theta,*rpy,**R,d,n; //use rpy[3]
  rpy=(double *)malloc(sizeof(double)*3);

  s=getCrossPro(v1,v2);
  getuv(s);
  theta=getAngleD(v1,v2)*pi/180;//acos(getDotPro(v1,v2)/(getMag(v1)*getMag(v2)));
  R=getR(s,theta);

  
  switch(ax){
    case 0:
      rpy[1]=asin(-1*R[2][0]);
      rpy[0]=asin(R[2][1]/cos(rpy[1]));
      rpy[2]=asin(R[1][0]/cos(rpy[1]));
      break;
    case 1:
      rpy[2]=asin(-1*R[0][1]);
      rpy[1]=asin(R[0][2]/cos(rpy[2]));
      rpy[0]=asin(R[2][1]/cos(rpy[2]));
      break;
    case 2:
      rpy[0]=asin(-1*R[1][2]);
      rpy[1]=asin(R[0][2]/cos(rpy[0]));
      rpy[2]=asin(R[1][0]/cos(rpy[0]));
      break;
  }
  

  for(i=0;i<3;i++)rpy[i]*=(180/pi);

  //mlfix
  free(s);
  for(i=0;i<3;i++)free(R[i]);
  free(R);
  //for(i=0;i<3;i++)printf("%f ",rpy[i]);printf("\n");
  return rpy;  ///degree
}

*/

double *getDm(double **cords,double *charge,int start_atom_no,int end_atom_no,double *point){
  int i,j;
  double *dm=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++)dm[i]=0;
  for(i=start_atom_no;i<=end_atom_no;i++){
    for(j=0;j<3;j++){
      //printf("%f: %f %f\n",point[j],cords[i][j],cords[i][j]-point[j]);
      dm[j]+=charge[i]*(cords[i][j]-point[j]);
    }
  }
  
  return dm;

}
