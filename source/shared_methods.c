#include "header_rotation.h"

double axis[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double *charge;

void init(char *file_name){
  int i=0;
  FILE *file = fopen(file_name,"r");
  char s[150];
  fgets(s,100,file);
  printf("atoms = %d\n",atoms=atoi(s));

  c = (char **)malloc(sizeof(char*)*atoms);

  for(i=0;i<atoms;i++){
    *(c+i) = (char *)malloc(sizeof(char)*2);
  }
  f1cords = (double **)malloc(sizeof(double*)*atoms);
  for(i=0;i<atoms;i++){
    *(f1cords+i) = (double *)malloc(sizeof(double)*3);
  }

  f2cords = (double **)malloc(sizeof(double*)*atoms);
  for(i=0;i<atoms;i++){
    *(f2cords+i) = (double *)malloc(sizeof(double)*3);
  }
 
  orig_f1cords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
    *(orig_f1cords+i)=(double *)malloc(sizeof(double)*3);
  }
  orig_f2cords=(double **)malloc(sizeof(double *)*atoms);
  for(i=0;i<atoms;i++){
    *(orig_f2cords+i)=(double *)malloc(sizeof(double)*3);
  }
 
 //charge
  charge = (double *)malloc(sizeof(double)*atoms);
  readCharge(charge_file_name);
  fclose(file);
}

void readCharge(char *file_name){
  FILE *file;
  int i,j;
  char s[150],*ch;

  file=fopen(file_name,"r");
  fgets(s,150,file);
  fgets(s,150,file);
  fgets(s,150,file);

  for(i=0;i<atoms&&fgets(s,150,file)!=NULL;i++){
  
    strtok(s," ");
    strtok(NULL," ");
    ch=strtok(NULL,"  ");
      //}else if(isdigit(ch[0])||ch[0]=='-'||ch[0]=='+'){
    charge[i]=atof(ch);
       
    //printf("%d) %s  %f\n",i,c[i],charge[i]);
  }
}

void printCords(double**cords,char**c){
  int i=0,j=0;
  for(i=0;i<atoms;i++){
    printf("%d) %s ",i,c[i]);
    for(j=0;j<3;j++){
      printf("%1.6f ",cords[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

int readFile(FILE *file,double** cords,char **c,int frame){
  int i=0,j=0,found=0;
  char s[150],*ch;
  fgets(s,150,file);
  while(frame>=0&&fgets(s,150,file)!=NULL){
    strtok(s," ");
    //strtok(NULL," ");
    if(frame==atoi(strtok(NULL," "))){
      found=1;
      break;
    }
    for(j=0;j<atoms+1&&fgets(s,150,file)!=NULL;j++);
    i++;
  }
  if(!found)return 0;

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
  return 1;
}

double getDist(double *p1,double *p2){
  return sqrt(pow(p1[0]-p2[0],2)+pow(p1[1]-p2[1],2)+pow(p1[2]-p2[2],2));
}

double getMag(double *a){
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

void getuv(double *v){
  int i;
  double d=getMag(v);
  if(d==0)return;
  for(i=0;i<3;i++){
    v[i]/=d;
  }
}

double getDotPro(double *a,double *b){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double *getCrossPro(double *a,double *b){
  double *crossPro=(double *)malloc(sizeof(double)*3);
  crossPro[0]=a[1]*b[2] - a[2]*b[1];
  crossPro[1]=-1*a[0]*b[2] + a[2]*b[0];
  crossPro[2]=a[0]*b[1] - a[1]*b[0];
  return crossPro;
}

double **getMatPro(double **a,int m1,int n1,double **b,int m2,int n2){
  int i,j,k;
  double **c;
  double d;
  if(n1!=m2){
    printf("Dimension Error\n");
    return NULL;
  }else{
    c=(double **)malloc(sizeof(double *)*m1);
    for(i=0;i<m1;i++){
      *(c+i)=(double *)malloc(sizeof(double)*n2);
    }
    for(i=0;i<m1;i++){
      for(j=0;j<n2;j++){d=0.0;
        for(k=0;k<n1;k++){
          d+=a[i][k]*b[k][j];
        }c[i][j]=d;
      }
    }
  return c;
  }
}

double *getCom(double **cords,int atoms_n){
  int i,j,k;
  double *p=(double*)malloc(sizeof(double)*3);
  double d,d1;
  for(i=0;i<3;i++){
    d=d1=0;
    for(j=0;j<atoms_n;j++){
      d+=cords[j][i]*getAtmass(j,c);
      d1+=getAtmass(j,c);
    }
    p[i]=d/d1;
  }
  return p;
}

double getAtmass(int n,char **c){
  if(c[n][0]=='C')return 12;
  else if(c[n][0]=='H')return 1;
  else if(c[n][0]=='N')return 14;
  else if(c[n][0]=='O')return 16;
  else if(c[n][0]=='S'&&c[n][1]=='i')return 28;
  else if(c[n][0]=='P')return 31;
  else if(c[n][0]=='F')return 18;
  else if(c[n][0]=='B')return 10;
  else printf("%d wtfffff\n",n);
  return -10000000000;
}

void shiftOrigin(double *origin,double *v){
  int i=0,j=0,k=0;
  double *s;
  double **R,**v1,**v2;
  //double axis[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  double d=0.0,d1=0.0,d2=0.0,theta=0.0;
  double align_v[3]={4,4,4},*atm_v;

  for(i=0;i<atoms;i++){
    for(j=0;j<3;j++){
      f1cords[i][j]-=origin[j];
      f2cords[i][j]-=origin[j];
    }
  }

  getuv(v);
  s=getCrossPro(v,axis[ax]);
  getuv(s);

  theta=acos(getDotPro(v,axis[ax]));
  R=getR(s,theta);
  for(i=0;i<atoms;i++){
    v1=(double**)malloc(sizeof(double*)*3);
    for(j=0;j<3;j++) *(v1+j)=(double*)malloc(sizeof(double));

    for(j=0;j<3;j++) v1[j][0]=f1cords[i][j];

    v2=getMatPro(R,3,3,v1,3,1);
    for(j=0;j<3;j++) f1cords[i][j]=v2[j][0];

    for(j=0;j<3;j++) v2[j][0]=f2cords[i][j];

    //mlfix
    for(k=0;k<3;k++)free(v1[k]);
    free(v1);

    v1=getMatPro(R,3,3,v2,3,1);
    for(j=0;j<3;j++) f2cords[i][j]=v1[j][0];
   
    //mlfix
    for(k=0;k<3;k++)free(v2[k]);
    free(v2);
    for(k=0;k<3;k++)free(v1[k]);
    free(v1);
  }

  //alignment: global atm is reference atom vector.....now changed to COM of whole system(not sure!)
  //atm_v=getCom(f2cords,207);//should be changed******
  //**atm_v = reference atom
  //**alig_v = reference direction
  atm_v=getCom(f2cords,atoms/2);//*
  //atm_v=(double *)malloc(sizeof(double)*3);for(i=0;i<3;i++)atm_v[i]=f1cords[0][i];
  align_v[ax]=atm_v[ax]=0;
  //for(i=0;i<3;i++)printf("%f %f\n",atm_v[i],align_v[i]);printf("\n");
  getuv(align_v);
  getuv(atm_v);

  //mlfix
  free(s);

  s=getCrossPro(atm_v,align_v);
  getuv(s);
  double sign=cos(acos(getDotPro(axis[ax],s)));
  theta=acos(getDotPro(align_v,atm_v));
  //printf("theta = %f sign = %f\n",theta*180/pi,sign);
  
  //mlfix
  for(i=0;i<3;i++)free(R[i]);
  free(R);

  R=getR(axis[ax],sign*theta);

  for(i=0;i<atoms;i++){
    v1=(double**)malloc(sizeof(double*)*3);
    for(j=0;j<3;j++) *(v1+j)=(double*)malloc(sizeof(double));

    for(j=0;j<3;j++) v1[j][0]=f1cords[i][j];

    v2=getMatPro(R,3,3,v1,3,1);
    for(j=0;j<3;j++) f1cords[i][j]=v2[j][0];

    for(j=0;j<3;j++) v2[j][0]=f2cords[i][j];

    //mlfix
    for(k=0;k<3;k++)free(v1[k]);
    free(v1);

    v1=getMatPro(R,3,3,v2,3,1);
    for(j=0;j<3;j++) f2cords[i][j]=v1[j][0];

    //mlfix
    for(k=0;k<3;k++)free(v2[k]);
    free(v2);
    for(k=0;k<3;k++)free(v1[k]);
    free(v1);
  }
  /*for(i=0;i<3;i++)atm_v[i]=f1cords[atm][i];
  atm_v[ax]=0;
  getuv(atm_v);
  d=acos(getDotPro(align_v,atm_v))*180/pi;
  //printf("%f %f %f %f\n\n",theta*180/pi,getDotPro(align_v,atm_v),d,(theta*180/pi)-d);
  */

  //mlfix
  free(s);
  for(i=0;i<3;i++)free(R[i]);
  free(R);
  //free(v1);
  free(atm_v);
  //free(axis);
}

double getAngleD(double *a,double *b){
  double d;
  d=getDotPro(a,b)/(getMag(a)*getMag(b));
  d=floor(d*1000000)/1000000;
  d=acos(d);
  return d*180/pi;
}

double getAngleR(double *a,double *b){
  double d;
  d=getDotPro(a,b)/(getMag(a)*getMag(b));
  d=floor(d*1000000)/1000000;
  d=acos(d);
  return d;
}

double dist(double **cords,int a,int b){
  return sqrt(pow(cords[a][0]-cords[b][0],2)+pow(cords[a][1]-cords[b][1],2)+pow(cords[a][2]-cords[b][2],2));
}

void writeOutput(FILE *file,int frame,double *angles){
  fprintf(file,"%d   %f\n",frame,angles[ax]);
}


//**************NO NEED TO CHANGE NAMES HERE*********//

void createFilez(FILE *f1,FILE *f2,int atoms,double **f1cords,double **f2cords,double p[],double p1[]){
  //FILE *f1;
  //FILE *f2;
  int i=0,j=0;
  /*char filename[50];
  sprintf(filename,"output/frame_%d.xyz",frame);
  f1 = fopen(filename,"w");
  sprintf(filename,"output/frame_%d.com",frame);
  f2 = fopen(filename,"w");
  */
  fprintf(f1,"%d\n\n",2*atoms+2);
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
  fprintf(f1,"B %f %f %f\n",p[0],p[1],p[2]);
  fprintf(f1,"B %f %f %f\n",p1[0],p1[1],p1[2]);
  //fprintf(f2,"\n");
  fprintf(f2,"%%nprocshared=24\n%%rwf=origi_o_2+.rwf\n%%nosave\n%%chk=origi_o_2+.chk\n# opt=maxcycles=500 rohf/genecp geom=connectivity\n\nMolden generated mol2\n\n0 1\n");
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f2,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f2,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f2,"%f ",f1cords[i][j]);
    }
    fprintf(f2,"\n");
  }
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f2,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f2,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f2,"%f ",f2cords[i][j]);
    }
    fprintf(f2,"\n");
  }
  fprintf(f2,"B %f %f %f\n",p[0],p[1],p[2]);
  fprintf(f2,"B %f %f %f\n",p1[0],p1[1],p1[2]);
  fprintf(f2,"\n");

  fclose(f1);
  fclose(f2);
}

void createFilezz(FILE *f1,FILE *f2,int atoms,int f1start,int f1end,int f2start,int f2end,double **f1cords,double **f2cords,double p[],double p1[]){
  if(f1start>=atoms||f1end>=atoms||f2start>=atoms||f2end>=atoms||f1end<f1start||f2end<f2start){
  printf("Check values createFilezz\n");
  return;
  }
  int i=0,j=0;
  fprintf(f1,"%d\n\n",(f1end-f1start+1)+(f2end-f2start+1)+2);
  for(i=f1start;i<=f1end;i++){
    if(c[i][1]!=0)
        fprintf(f1,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f1,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f1,"%f ",f1cords[i][j]);
    }
    fprintf(f1,"\n");
  }
  for(i=f2start;i<=f2end;i++){
    if(c[i][1]!=0)
        fprintf(f1,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f1,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f1,"%f ",f2cords[i][j]);
    }
    fprintf(f1,"\n");
  }
  fprintf(f1,"B %f %f %f\n",p[0],p[1],p[2]);
  fprintf(f1,"B %f %f %f\n",p1[0],p1[1],p1[2]);

  fprintf(f2,"%%nprocshared=24\n%%rwf=origi_o_2+.rwf\n%%nosave\n%%chk=origi_o_2+.chk\n# opt=maxcycles=500 rohf/genecp geom=connectivity\n\nMolden generated mol2\n\n0 1\n");
  for(i=f1start;i<=f1end;i++){
    if(c[i][1]!=0)
        fprintf(f2,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f2,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f2,"%f ",f1cords[i][j]);
    }
    fprintf(f2,"\n");
  }
  for(i=f2start;i<=f2end;i++){
    if(c[i][1]!=0)
        fprintf(f2,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f2,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f2,"%f ",f2cords[i][j]);
    }
    fprintf(f2,"\n");
  }
  fprintf(f2,"B %f %f %f\n",p[0],p[1],p[2]);
  fprintf(f2,"B %f %f %f\n",p1[0],p1[1],p1[2]);
  fprintf(f2,"\n");

  fclose(f1);
  fclose(f2);
}

//***********CHANGE THE NAME OF OUTPUT FILES DOWN HERE*****************//
/*
void createFile(char **c,int atoms,double theta,double **rcords){
  FILE *f1;
  FILE *f2;
  int i=0,j=0;

  char filename[100];
  if(theta<0){
    sprintf(filename,"result_rotz_neg_%.0f.xyz",-1*theta);
  }else{
    sprintf(filename,"result_rotz_pos_%.0f.xyz",theta);
  }
  f1 = fopen(filename,"w");

  if(theta<0){
    sprintf(filename,"result_rotz_neg_%.0f.com",-1*theta);
  }else{
    sprintf(filename,"result_rotz_pos_%.0f.com",theta);
  }
  f2 = fopen(filename,"w");

  fprintf(f1,"%d\n\n",atoms);
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f1," %c%c ",c[i][0],c[i][1]);
    else
        fprintf(f1," %c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f1,"%f ",rcords[i][j]);
    }
    fprintf(f1,"\n");
  }
  fprintf(f1,"\n");

  fprintf(f2,"%%nprocshared=24\n%%rwf=origi_o_2+.rwf\n%%nosave\n%%chk=origi_o_2+.chk\n# opt=maxcycles=500 rohf/genecp geom=connectivity\n\nMolden generated mol2\n\n0 1\n");
  for(i=0;i<atoms;i++){
    if(c[i][1]!=0)
        fprintf(f2,"%c%c ",c[i][0],c[i][1]);
    else
        fprintf(f2,"%c%c ",c[i][0],' ');
    for(j=0;j<3;j++){
      fprintf(f2,"%f ",rcords[i][j]);
    }
    fprintf(f2,"\n");
  }
  fprintf(f2,"\n");

  fclose(f1);
  fclose(f2);
}
*/

