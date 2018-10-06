#include "source/header_rotation.h"

int atoms;
double **f1cords,**f2cords,**orig_f1cords,**orig_f2cords;
char **c;

int main(){
  int i=0,j=0;
  char *files[]={file1_name,file2_name};
  double *net_rpy,*rpy;
  init(files[0]);

  //********DIRECTIONALITY********//
  if(0){
  //step_size=getBestStepSize(files,1,sqrt(end_frame));
  printf("Best Step Size = %d\n",step_size);
  net_rpy=getDirectionality(files,start_frame,end_frame,step_size);
  if(net_rpy!=NULL){
    printf("\naxis = %d\n",ax);
    for(i=0;i<3;i++)printf("%f\n",net_rpy[i]);printf("\n");
  }else{
    printf("NULL angles\n");
  }
  //*/
  }
  //*********Direct Rotation b/w start and end frame***********//
  if(0){
  //ax=0;
  rpy=getRingRotation(files,start_frame,end_frame,1);
  if(rpy!=NULL){
    for(i=0;i<3;i++)printf("%f\n",rpy[i]);printf("\n");
  }else{
    printf("NULL angles\n");
  }
 
  //memory_leak_fix
  free(net_rpy);
  free(rpy);
  }
 

  //********CODE VERIFICATION***********//
 
  if(0){
  //********old rpy test***********//
  //testRPY(files,0,999999); 
 
  //********new rpy test**********//
  char *test_file_name="input/zz.xyz";
  char *test2files[]={test_file_name};
  init(test2files[0]);
  rpyTest2(test_file_name,files);
  //testRPY(files,0,1);
  }


//*******DIPOLE MOMENT*********//
  if(1){
  int f1=0,isf1;
  double *dm;
  //for(i=0;i<atoms;i++)printf("%d) %f\n",i,charge[i]); 
  FILE *file1=fopen(f1<=file_limit?files[0]:files[1],"r");
  isf1=readFile(file1,f1cords,c,f1);
  if(!isf1){
    printf("Error in reading file\n");
    return NULL;
  }else{
    printf("Found frame %d\n",f1);
  }


  double origin[3]={0,0,0};
  double *com=getCom(f1cords,atoms);
  //for(i=0;i<3;i++)printf("-->%f ",com[i]);printf("\n");
  dm=getDm(f1cords,charge,0,atoms-1,origin);//getCom(f1cords,atoms));
  //for(i=0;i<3;i++)printf("%f ",dm[i]);printf("\n");
  printf("%f\n",getMag(dm));
  
  FILE *dm_file=fopen("output/dm.csv","w");
  fprintf(dm_file,"charges,x,y,z\n");
  for(i=0;i<atoms;i++){
    fprintf(dm_file,"%f,%f,%f,%f\n",charge[i],f1cords[i][0],f1cords[i][1],f1cords[i][2]);
  }
  
  fclose(dm_file);
  //dm ring rotation
  dmGetRingRotation(files,0,5000,1);
  getRingRotation(files,0,5000,1);
  /*net_rpy=getDirectionality(files,0,10000,step_size);
  if(net_rpy!=NULL){
    printf("\naxis = %d\n",ax);
    for(i=0;i<3;i++)printf("%f\n",net_rpy[i]);printf("\n");
  }else{
    printf("NULL angles\n");
  }

  free(net_rpy);*/
  net_rpy=dmGetDirectionality(files,0,10000,step_size);
  if(net_rpy!=NULL){
    printf("\naxis = %d\n",ax);
    for(i=0;i<3;i++)printf("%f\n",net_rpy[i]);printf("\n");
  }else{
    printf("NULL angles\n");
  }

  }
  

  //************DEBUG*************//
      
  if(0){
  //ax=0;
  int f1=0,f2=8711;
  double *f1com,*f2com;
  rpy=getRingRotation(files,f1,f2,1);
  if(rpy!=NULL){
    for(i=0;i<3;i++)printf("%f\n",rpy[i]);printf("\n");
  }else{
    printf("NULL angles\n");
  }
  
  checkTransformation();
  
  char filename[50];
  sprintf(filename,"output/shfited_%d-%d.xyz",f1,f2);
  FILE *file1 = fopen(filename,"w");
  sprintf(filename,"output/shifted_%d-%d.com",f1,f2);
  FILE *file2 = fopen(filename,"w");
  translateMol(axis[ax],5,0,atoms-1,f2cords);
  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(f2cords,ring_atom);
  //createFilez(file1,file2,atoms,f1cords,f2cords,f1com,f2com);
  createFilezz(file1,file2,atoms,0,ring_atom-1,1,10,f1cords,f2cords,f1com,f2com);
  free(f1com);
  free(f2com);

  sprintf(filename,"output/orig_%d-%d.xyz",f1,f2);
  file1 = fopen(filename,"w");
  sprintf(filename,"output/orig_%d-%d.com",f1,f2);
  file2 = fopen(filename,"w");
  f1com=getCom(orig_f1cords,ring_atom);
  f2com=getCom(orig_f2cords,ring_atom);
  createFilez(file1,file2,atoms,orig_f1cords,orig_f2cords,f1com,f2com);

  free(f1com);
  free(f2com);
  sprintf(filename,"output/frames_%d.xyz",f1);
  file1 = fopen(filename,"w");
  sprintf(filename,"output/frames_%d.com",f1);
  file2 = fopen(filename,"w");
  f1com=getCom(f1cords,ring_atom);
  f2com=getCom(orig_f1cords,ring_atom);
  createFilez(file1,file2,atoms,orig_f1cords,f1cords,f1com,f2com);

  free(f1com);
  free(f2com);
  sprintf(filename,"output/frames_%d.xyz",f2);
  file1 = fopen(filename,"w");
  sprintf(filename,"output/frames_%d.com",f2);
  file2 = fopen(filename,"w");
  f1com=getCom(f2cords,ring_atom);
  f2com=getCom(orig_f2cords,ring_atom);
  createFilez(file1,file2,atoms,orig_f2cords,f2cords,f1com,f2com);

  

  //memory_leak_fix
  free(rpy);
  free(f1com);
  free(f2com);
  for(i=0;i<atoms;i++)free(orig_f1cords[i]);
  free(orig_f1cords);
  for(i=0;i<atoms;i++)free(orig_f2cords[i]);
  free(orig_f2cords);
  }
  
  
  return 0;
}



