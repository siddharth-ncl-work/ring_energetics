#ifndef ROTATION_H
#define ROTATION_H

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>

#define pi acos(-1)

void init(char*);
int readFile(FILE *,double **,char **,int);
void printCords(double **,char **);
double * getRingRotation(char *[],int,int,int);
double *getDirectionality(char *[],int,int,int);
void shiftOrigin(double *,double *);
double *getCom(double **,int);
double getAtmass(int,char **);
double getDist(double *,double *);
double **getMatPro(double **,int,int,double **,int,int);
void getuv(double *);
double getDotPro(double *,double *);
double getMag(double *);
double **getR(double *,double);
double *getRPY(double *,double *);
double getAngleD(double *,double *);
double dist(double **,int,int);
//void createFile(char **,int,double,double **);
void writeOutput(FILE *,int,double *);
double *getCrossPro(double *,double *);
void createFilez(FILE *,FILE *,int,double **,double **,double [],double []);
void createFilezz(FILE *,FILE *,int,int,int,int,int,double **,double **,double [],double []);
int getBestStepSize(char *[],int,int);
double * getTrackRotation(char *[],int,int,int,int,int); //track rotation

//*****Torque****//
void readCharge(char *);
double *getTorqueInFrame(int);
double *getElectroForce(int);
void writeTorqueOutput(FILE *,int,double *,double);
//******translation********//
double getNetTranslation(char *[],int,int,int);
double getRingTranslation(char *[],int,int,int);
double *getTransformationTrans(double [],double []);
double getTrans(double [],double []);
//****test******//
void mol_rotx(double **,int,double,double**);
void mol_roty(double **,int,double,double**);
void mol_rotz(double **,int,double,double**);
void rot(double *,double *,double *);
void rot1(double *,double *,double *,double *);
double *getuv1(double *);
void testRPY(char *[],int,int);
//IMP test2
void test2ReadFile(FILE *,double **,char **);
void test2Rot(double *,double);
void test2Trans(double *,double);
void createTestFrames();
double *getAxis();
void rpyTest2withTrans();
void avgRpyTest2withTrans(double *,double *);
void rpyTest2(char *,char *[]);
//debug
double *translate(double *,double,double *);
void translateMol(double *,double,int,int,double **);
void checkTransformation();
double frameAtomDist(double **,double **,int);
//****dipole_moment****//
double *dmGetDirectionality(char *[],int,int,int);
double *dmGetRingRotation(char **,int,int,int);
double *dmGetTrackRotation(char **,int,int,int,int,int);
double *getDm(double **,double *,int,int,double *);


//gobal variables
extern int start_frame;
extern int end_frame;
extern int  step_size;           //****the old "inc" variable
extern char *file1_name;
extern char *file2_name;
extern int file_limit;
extern int ring_atom;  //number of ring atoms
extern int track_start_atom_no;
extern int track_end_atom_no;
extern int atoms,ax;
extern double **f1cords,**f2cords,**orig_f1cords,**orig_f2cords;
extern char **c;
//torque
extern double axis[3][3];
extern double *charge;
extern char *charge_file_name;
extern double K;
#endif
