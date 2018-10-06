#include "header_rotation.h"


double *getTorqueInFrame(int frame_no){
  int i,j;
  double r[3],*f,*torque,*net_torque,t,*ring_com;
  net_torque=(double *)malloc(sizeof(double)*3);
  ring_com=getCom(f2cords,ring_atom);
  for(i=0;i<ring_atom;i++){
    f=getElectroForce(i);
    for(j=0;j<3;j++){
      r[j]=f2cords[i][j]-ring_com[j];
    }
    torque=getCrossPro(r,f);
    for(j=0;j<3;j++){
      net_torque[j]+=torque[j];
    }
  }
  //t=getAngleD(axis[ax],net_torque);
  //printf("theta %f\n",t); 
  free(f);
  return net_torque;
}

void writeTorqueOutput(FILE *file,int frame,double *net_torque,double angle){
  int tt,dd;
  double t=getAngleD(axis[ax],net_torque);
  printf("torque = %f\n",getMag(net_torque));
  printf("theta %f\n",t);
  if(t>90){
    tt=-1;
  }else{
    tt=1;
  }
  if(angle<0){
    dd=-1;
  }else{
    dd=1;
  }
  fprintf(file,"%d   %d\n",tt,dd);
  free(net_torque);
}

double *getElectroForce(int target){
  int i,j;
  double *netf=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++){
    for(j=ring_atom;j<atoms;j++){
      netf[i]+=charge[j]*(f2cords[target][i]-f2cords[j][i])/pow(dist(f2cords,j,target),3);
    }
    netf[i]*=charge[target]*K;
  }
  return netf;  
}
