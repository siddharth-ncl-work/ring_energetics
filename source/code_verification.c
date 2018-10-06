#include "header_rotation.h"

//old rpy test
void testRPY(char *files[],int f1,int f2){
  int i=0,j=0,isf1=0,isf2=0;
  double rpy[3]={0,0,0},*v,**rcords,v1[3]={10,0,-33.5};//*t;
  FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
  FILE *file2=fopen(f2<=file_limit?files[0]:files[1],"r");
  rcords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++) *(rcords+i)=(double *)malloc(sizeof(double)*3);

  isf1=readFile(file1,f1cords,c,f1);
  isf2=readFile(file2,f2cords,c,f2);
  if(!isf1||!isf2){
    printf("Error in reading files\n");
    return;
  }else{
    printf("Found frame %d and %d\n",f1,f2);
  }

//#single vector

  double t[]={0,0,0};
  double *v2=(double *)malloc(sizeof(double)*3);
  rot(v1,&t[0],v2);
  v=getRPY(v1,v2);
  double *v3=(double *)malloc(sizeof(double)*3);
  rot(v1,v,v3);
  for(i=0;i<3;i++){
    printf("%f %f %f\n",v[i],v2[i],v3[i]);
  }


//#Molecule
/*
  mol_rotx(f1cords,atoms,90,rcords);
  mol_roty(rcords,atoms,180,f2cords);
  mol_rotz(f2cords,atoms,270,rcords);

  for(i=0;i<atoms;i++){
    v=getRPY(f1cords[i],rcords[i]);
    for(j=0;j<3;j++)printf("%f ",v[j]);printf("\n"); *
    rot(f1cords[i],v,f2cords[i]);
    for(j=0;j<3;j++){
      rpy[j]+=v[j];
    }
  }
*/


//#frames
/*
  for(i=0;i<atoms;i++){///degrees
    v=getRPY(f1cords[i],f2cords[i]);
    printf("%d)\n",i);
    for(j=0;j<3;j++)printf("%f ",f1cords[i][j]); printf("\n");
    for(j=0;j<3;j++)printf("%f ",f2cords[i][j]); printf("\n");
    for(j=0;j<3;j++)printf("%f ",v[j]); printf("\n");
    t=getTrans(f1cords[i],f2cords[i]);
    //rot(f1cords[i],v,rcords[i]);
    rot1(f1cords[i],v,t,rcords[i]);
    for(j=0;j<3;j++){
      rpy[j]+=v[j];
    }
  }
*/

//#Average
  /*
  printf("RPY = ");
  for(i=0;i<3;i++){
    printf("%f ",rpy[i]);
    printf("%f, ",rpy[i]/=atoms);
  }printf("\n\n");

  mol_rotx(f1cords,atoms,rpy[0],f2cords);
  mol_roty(f2cords,atoms,rpy[1],f2cords);
  mol_rotz(f2cords,atoms,rpy[2],f2cords);
  */

//#Don't know where does this belongs
/*
  for(i=0;i<atoms;i++){
    printf("%d)  ",i);
    for(j=0;j<3;j++){
      printf("%f -->%f %f   ",f1cords[i][j],rcords[i][j],f2cords[i][j]);
    }printf("\n");
  }
*/

  fclose(file1);
  fclose(file2);
}



void mol_rotx(double **cords,int atoms,double theta1,double **rcords){

  int i = 0, j = 0;

  double theta = 0.0;

  double rotm[3][3];
    theta = theta1*acos(-1)/180;
    rotm[0][0] = 1;
    rotm[0][1] = 0;
    rotm[0][2] = 0;
    rotm[1][0] = 0;
    rotm[1][1] = cos(theta);
    rotm[1][2] = -1*sin(theta);
    rotm[2][0] = 0;
    rotm[2][1] = sin(theta);
    rotm[2][2] = cos(theta);


  for(i=0;i<atoms;i++){
    for(j=0;j<3;j++){
      rcords[i][j] = rotm[j][0]*cords[i][0]+rotm[j][1]*cords[i][1]+rotm[j][2]*cords[i][2];
    }
  }

}

void mol_roty(double **cords,int atoms,double theta1,double **rcords){

  int i = 0, j = 0;

  double theta = 0.0;

  double rotm[3][3];
    theta = theta1*acos(-1)/180;
    rotm[0][0] = cos(theta);
    rotm[0][1] = 0;
    rotm[0][2] = sin(theta);
    rotm[1][0] = 0;
    rotm[1][1] = 1;
    rotm[1][2] = 0;
    rotm[2][0] = -1*sin(theta);
    rotm[2][1] = 0;
    rotm[2][2] = cos(theta);


  for(i=0;i<atoms;i++){
    for(j=0;j<3;j++){
      rcords[i][j] = rotm[j][0]*cords[i][0]+rotm[j][1]*cords[i][1]+rotm[j][2]*cords[i][2];
    }
  }

}

void mol_rotz(double **cords,int atoms,double theta1,double **rcords){

  int i = 0, j = 0;

  double theta = 0.0;

  double rotm[3][3];
    theta = theta1*acos(-1)/180;
    rotm[0][0] = cos(theta);
    rotm[0][1] = -1*sin(theta);
    rotm[0][2] = 0;
    rotm[1][0] = sin(theta);
    rotm[1][1] = cos(theta);
    rotm[1][2] = 0;
    rotm[2][0] = 0;
    rotm[2][1] = 0;
    rotm[2][2] = 1;


  for(i=0;i<atoms;i++){
    for(j=0;j<3;j++){
      rcords[i][j] = rotm[j][0]*cords[i][0]+rotm[j][1]*cords[i][1]+rotm[j][2]*cords[i][2];
    }
  }

}

void rot(double *v1,double *t,double *v2){
  int i=0,j=0;
  double rad=pi/180;
  double tx=t[0]*rad,ty=t[1]*rad,tz=t[2]*rad,**c,**Rx,**Ry,**Rz;

  for(i=0;i<3;i++)printf("%f ",t[i]);printf("\n");
  Rx=(double **)malloc(sizeof(double *)*3);
  Ry=(double **)malloc(sizeof(double *)*3);
  Rz=(double **)malloc(sizeof(double *)*3);
  for(i=0;i<3;i++){
    *(Rx+i)=(double *)malloc(sizeof(double)*3);
    *(Ry+i)=(double *)malloc(sizeof(double)*3);
    *(Rz+i)=(double *)malloc(sizeof(double)*3);
  }

  Rx[0][0] = 1;
  Rx[0][1] = 0;
  Rx[0][2] = 0;
  Rx[1][0] = 0;
  Rx[1][1] = cos(tx);
  Rx[1][2] = -1*sin(tx);
  Rx[2][0] = 0;
  Rx[2][1] = sin(tx);
  Rx[2][2] = cos(tx);

  Ry[0][0] = cos(ty);
  Ry[0][1] = 0;
  Ry[0][2] = sin(ty);
  Ry[1][0] = 0;
  Ry[1][1] = 1;
  Ry[1][2] = 0;
  Ry[2][0] = -1*sin(ty);
  Ry[2][1] = 0;
  Ry[2][2] = cos(ty);

  Rz[0][0] = cos(tz);
  Rz[0][1] = -1*sin(tz);
  Rz[0][2] = 0;
  Rz[1][0] = sin(tz);
  Rz[1][1] = cos(tz);
  Rz[1][2] = 0;
  Rz[2][0] = 0;
  Rz[2][1] = 0;
  Rz[2][2] = 1;

  c=(double **)malloc(sizeof(double *)*3);
  for(i=0;i<3;i++){
    *(c+i)=(double *)malloc(sizeof(double));
    c[i][0]=v1[i];
  }
  c=getMatPro(Rx,3,3,c,3,1);
  c=getMatPro(Ry,3,3,c,3,1);
  c=getMatPro(Rz,3,3,c,3,1);

  for(i=0;i<3;i++){
    v2[i]=c[i][0];
  }
}

void rot1(double *v1,double *t,double *trans,double *v2){
  int i=0,j=0;
  double rad=pi/180;
  double tx=t[0]*rad,ty=t[1]*rad,tz=t[2]*rad,**c,**Rx,**Ry,**Rz,**Rt;

  Rx=(double **)malloc(sizeof(double *)*3);
  Ry=(double **)malloc(sizeof(double *)*3);
  Rz=(double **)malloc(sizeof(double *)*3);
  Rt=(double **)malloc(sizeof(double *)*4);
  for(i=0;i<3;i++){
    *(Rx+i)=(double *)malloc(sizeof(double)*3);
    *(Ry+i)=(double *)malloc(sizeof(double)*3);
    *(Rz+i)=(double *)malloc(sizeof(double)*3);
    *(Rt+i)=(double *)malloc(sizeof(double)*4);
  }
  *(Rt+i)=(double *)malloc(sizeof(double)*4);

  Rx[0][0] = 1;
  Rx[0][1] = 0;
  Rx[0][2] = 0;
  Rx[1][0] = 0;
  Rx[1][1] = cos(tx);
  Rx[1][2] = -1*sin(tx);
  Rx[2][0] = 0;
  Rx[2][1] = sin(tx);
  Rx[2][2] = cos(tx);

  Ry[0][0] = cos(ty);
  Ry[0][1] = 0;
  Ry[0][2] = sin(ty);
  Ry[1][0] = 0;
  Ry[1][1] = 1;
  Ry[1][2] = 0;
  Ry[2][0] = -1*sin(ty);
  Ry[2][1] = 0;
  Ry[2][2] = cos(ty);

  Rz[0][0] = cos(tz);
  Rz[0][1] = -1*sin(tz);
  Rz[0][2] = 0;
  Rz[1][0] = sin(tz);
  Rz[1][1] = cos(tz);
  Rz[1][2] = 0;
  Rz[2][0] = 0;
  Rz[2][1] = 0;
  Rz[2][2] = 1;

  Rz=getMatPro(Rz,3,3,Ry,3,3);
  Rz=getMatPro(Rz,3,3,Rx,3,3);
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      Rt[i][j]=Rz[i][j];
    }
  }
  for(i=0;i<3;i++) Rt[i][3]=trans[i];
  for(i=0;i<3;i++) Rt[3][i]=0;
  Rt[3][3]=1;

  /*
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      printf("%f ",Rt[i][j]);
    }printf("\n");
  }printf("\n");
  */

  c=(double **)malloc(sizeof(double *)*4);
  for(i=0;i<3;i++){
    *(c+i)=(double *)malloc(sizeof(double));
    c[i][0]=v1[i];
  }
  *(c+i)=(double *)malloc(sizeof(double));
  c[i][0]=1;

  c=getMatPro(Rt,4,4,c,4,1);
  for(i=0;i<3;i++){
    v2[i]=c[i][0];
  }
}
                                                            
double *getuv1(double *v){
  int i=0;
  double d=getMag(v);
  double *uv=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++){
    uv[i]=v[i]/d;
  }
  return uv;
}

/*
double *getTrans(double v1[],double v2[]){
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
*/

//*****__new rpy test__IMP****//

void test2ReadFile(FILE *file,double **cords,char** c){
int i=0,j=0;//,found=0;
  char s[150],*ch;
  fgets(s,150,file);
  fgets(s,150,file);
  /*
  while(frame>=0&&fgets(s,150,file)!=NULL){
    strtok(s," ");
    if(frame==atoi(strtok(NULL," "))){
      found=1;
      break;
    }
    for(j=0;j<atoms+1&&fgets(s,150,file)!=NULL;j++);
    i++;
  }
  if(!found)return 0;
  */
  for(i=0;i<atoms&&fgets(s,150,file)!=NULL;i++){
    j=0;
    ch=strtok(s," ");
    while(ch!=NULL){
      if(isalpha(ch[0])){
        strcpy(*(c+i),ch);
      }else if(isdigit(ch[0])||ch[0]=='-'||ch[0]=='+'){
         cords[i][j++]=atof(ch);
      }
      ch=strtok(NULL,"  ");
    }
    //printf("%d) %s  %f %f %f\n",i,c[i],cords[i][0],cords[i][1],cords[i][2]);
  }
  //return 1;  
}

void test2Rot(double *s,double t){
  int i,j,k;
  double **R,**v1,**v2;
  R=getR(s,t);
  for(i=0;i<ring_atom;i++){
    v1=(double**)malloc(sizeof(double*)*3);
    for(j=0;j<3;j++) *(v1+j)=(double*)malloc(sizeof(double));

    for(j=0;j<3;j++) v1[j][0]=f1cords[i][j];

    v2=getMatPro(R,3,3,v1,3,1);
    for(j=0;j<3;j++) f2cords[i][j]=v2[j][0];
     
    //for(k=0;k<3;k++)printf("%f %f\n",f1cords[i][k],f2cords[i][k]);
    printf("-->%f\n",getAngleD(f1cords[i],f2cords[i]));

    //mlfix
    for(k=0;k<3;k++)free(v2[k]);
    free(v2);
    for(k=0;k<3;k++)free(v1[k]);
    free(v1);
  }

}

void test2Trans(double *s,double d){
  int i,j;
  for(i=0;i<ring_atom;i++){
    for(j=0;j<3;j++){
      f2cords[i][j]+=d*s[j];
    }
  }

}

void createTestFrames(){
  FILE *f1;
  int i,j;
  f1=fopen("input/test_frames.xyz","w");
  fprintf(f1,"%d\n",atoms);
  fprintf(f1,"frame 0 xyzgerte \n");
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f1,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f1,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f1,"%f ",f1cords[i][j]);
    }
    fprintf(f1,"\n");
  }

  fprintf(f1,"%d\n",atoms);
  fprintf(f1,"frame 1 xyszdf\n");
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f1,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f1,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f1,"%f ",f2cords[i][j]);
    }
    fprintf(f1,"\n");
  }
  fclose(f1);
}

double *getAxis(){
  int i=0;
  double *f1com,*f2com;
  double *ax=(double *)malloc(sizeof(double)*3);
  ax[0]=1;
  ax[1]=ax[2]=0;
  //f1com=getCom(f1cords,ring_atom);
  //f2com=getCom(f2cords,ring_atom);
  //for(i=0;i<3;i++)ax[i]=f2com[i]-f1com[i];
   
  getuv(ax);
  return ax;
} 

void rpyTest2withTrans(){
  int i,j;
  double *v,*t,**rcords;
  rcords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++) *(rcords+i)=(double *)malloc(sizeof(double)*3);
  /*
  for(i=0;i<atoms;i++){///degrees
    v=getRPY(f1cords[i],f2cords[i]);
    t=getTrans(f1cords[i],f2cords[i]);
    rot1(f1cords[i],v,t,rcords[i]);

    //mlfix
    free(v);
    free(t);
  }
  
  for(i=0;i<atoms;i++){
    printf("%d)  ",i);
    for(j=0;j<3;j++){
      printf("%f -->%f %f   ",f1cords[i][j],rcords[i][j],f2cords[i][j]);
    }printf("\n");
  }
  */
}

void avgRpyTest2withTrans(double *rpy,double *trans){
  int i=0,j=0,k=0;
  double rad=pi/180;
  double tx=rpy[0]*rad,ty=rpy[1]*rad,tz=rpy[2]*rad,**c,**Rx,**Ry,**Rz,**Rt,**v2,**rcords;

  rcords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
    *(rcords+i)=(double *)malloc(sizeof(double)*3);
    for(j=0;j<3;j++)rcords[i][j]=f1cords[i][j];
  }
  
  
  Rx=(double **)malloc(sizeof(double *)*3);
  Ry=(double **)malloc(sizeof(double *)*3);
  Rz=(double **)malloc(sizeof(double *)*3);
  Rt=(double **)malloc(sizeof(double *)*4);
  for(i=0;i<3;i++){
    *(Rx+i)=(double *)malloc(sizeof(double)*3);
    *(Ry+i)=(double *)malloc(sizeof(double)*3);
    *(Rz+i)=(double *)malloc(sizeof(double)*3);
    *(Rt+i)=(double *)malloc(sizeof(double)*4);
  }
  *(Rt+i)=(double *)malloc(sizeof(double)*4);

  Rx[0][0] = 1;
  Rx[0][1] = 0;
  Rx[0][2] = 0;
  Rx[1][0] = 0;
  Rx[1][1] = cos(tx);
  Rx[1][2] = -1*sin(tx);
  Rx[2][0] = 0;
  Rx[2][1] = sin(tx);
  Rx[2][2] = cos(tx);

  Ry[0][0] = cos(ty);
  Ry[0][1] = 0;
  Ry[0][2] = sin(ty);
  Ry[1][0] = 0;
  Ry[1][1] = 1;
  Ry[1][2] = 0;
  Ry[2][0] = -1*sin(ty);
  Ry[2][1] = 0;
  Ry[2][2] = cos(ty);

  Rz[0][0] = cos(tz);
  Rz[0][1] = -1*sin(tz);
  Rz[0][2] = 0;
  Rz[1][0] = sin(tz);
  Rz[1][1] = cos(tz);
  Rz[1][2] = 0;
  Rz[2][0] = 0;
  Rz[2][1] = 0;
  Rz[2][2] = 1;

  Rz=getMatPro(Rz,3,3,Ry,3,3);
  Rz=getMatPro(Rz,3,3,Rx,3,3);
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      Rt[i][j]=Rz[i][j];
    }
  }
  for(i=0;i<3;i++) Rt[i][3]=trans[i];
  for(i=0;i<3;i++) Rt[3][i]=0;
  Rt[3][3]=1;

  /*
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      printf("%f ",Rt[i][j]);
    }printf("\n");
 }printf("\n");
 */


  for(k=0;k<ring_atom;k++){
    c=(double **)malloc(sizeof(double *)*4);
    for(i=0;i<3;i++){
      *(c+i)=(double *)malloc(sizeof(double));
      c[i][0]=f1cords[k][i];
    }
    *(c+i)=(double *)malloc(sizeof(double));
    c[i][0]=1;

    v2=getMatPro(Rt,4,4,c,4,1);
    for(i=0;i<3;i++){
      rcords[k][i]=v2[i][0];
    }
   
    //mlfix
    for(i=0;i<4;i++)free(c[i]);
    free(c);
    for(i=0;i<4;i++)free(v2[i]);
    free(v2);
  }

  for(i=0;i<atoms;i++){
    printf("%d)  ",i);
    for(j=0;j<3;j++){
      printf("%f -->%f %f   ",f1cords[i][j],rcords[i][j],f2cords[i][j]);
    }printf("\n");
  }

  for(i=0;i<atoms;i++)free(rcords[i]);
  free(rcords); 
}

void rpyTest2(char *test_mol_file_name,char *files[]){
  int i,j;
  double *s,*f1com,*rpy;
  s=(double *)malloc(sizeof(double)*3);
  FILE *test_file=fopen(test_mol_file_name,"r");
  test2ReadFile(test_file,f1cords,c);
  //printCords(f1cords,c);
  for(i=0;i<atoms;i++){
    for(j=0;j<3;j++){
      f2cords[i][j]=f1cords[i][j];
    }
  }

  f1com=getCom(f1cords,ring_atom);
  for(i=0;i<3;i++){
    s[i]=f1cords[51][i]-f1com[i];//y.xyz-->40,24  zz.xyz-->30,52
  } 
  getuv(s);   
  //printf("-->%f\n",getMag(s));

  //f1com=getCom(f1cords,ring_atom);

  shiftOrigin(f1com,s);

  double *ax=getAxis();

  //test2Trans(ax,1);
  
  f1com=getCom(f1cords,ring_atom);
  for(i=0;i<3;i++){
    s[i]=f1cords[51][i]-f1com[i];//y.xyz-->40,24  zz.xyz-->30,52
  }
  getuv(s);
  printf("-->%f\n",getAngleD(ax,s)); 
  


  test2Rot(ax,60*pi/180);
 
  
  printf("-->>%f\n",getAngleD(f1cords[0],f2cords[0]));
  f1com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++){
    s[i]=f1cords[51][i]-f1com[i];//y.xyz-->40,24  zz.xyz-->30,52
  }
  getuv(s);
  printf("-->%f\n",getAngleD(ax,s));



  test2Trans(ax,5);


  f1com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++){
    s[i]=f1cords[51][i]-f1com[i];//y.xyz-->40,24  zz.xyz-->30,52
  }
  getuv(s);
  printf("-->%f\n",getAngleD(ax,s));

  f1com=getCom(f1cords,ring_atom);
  double *f2com=getCom(f2cords,ring_atom);
  for(i=0;i<3;i++){
    s[i]=f2com[i]-f1com[i];//y.xyz-->40,24  zz.xyz-->30,52
  }
  getuv(s);
  printf("-->>>%f\n",getAngleD(ax,s));

  printf("-->%f\n",getAngleD(f1cords[0],f2cords[0]));
  
  


  createTestFrames();
  //*rpyTest2withTrans();
  rpy=getRingRotation(files,0,1,1);
  //avgRpyTest2withTrans(rpy,trans);
  for(i=0;i<3;i++)printf("rpy %f\n",rpy[i]);
} 


//****DEBUG****//
//visualization stuff

double *translate(double *s,double n,double *v){
  int i,j;
  double *t=(double *)malloc(sizeof(double)*3);
  for(i=0;i<3;i++)t[i]=0;
  for(i=0;i<3;i++)t[i]=v[i]+n*s[i];
  return t;
}

void translateMol(double *s,double n,int start_atom_no,int end_atom_no,double **cords){
  int i,j;
  double *t;

  for(i=start_atom_no;i<=end_atom_no&&i<atoms;i++){
    t=translate(s,n,cords[i]);
    for(j=0;j<3;j++)cords[i][j]=t[j];
    free(t);
  }

}

//******CHECK TRANSFORMATION******//
void checkTransformation(){
  int i,j;
  double d,d1;
  for(i=0;i<atoms;i++){
    d=frameAtomDist(orig_f1cords,orig_f2cords,i);
    d=roundf(d * pow(10,10)) / pow(10,10);
    d1=frameAtomDist(f1cords,f2cords,i);
    d1=roundf(d1 * pow(10,10)) / pow(10,10);
    //printf("%d) %f %f\n",i,d,d1);
    if(d==d1){
      
    }else{
      printf("oops! %d\n",i);
    }
  }
}

double frameAtomDist(double **f1cords,double **f2cords,int n){
  return sqrt(pow(f1cords[n][0]-f2cords[n][0],2)+pow(f1cords[n][1]-f2cords[n][1],2)+pow(f1cords[n][2]-f2cords[n][2],2));
}

