//** file for all global variables **//

int ax=1;
char *file1_name="/home/vanka/siddharth/ruchi/code_file_v1.0.5_with_counter_ions_8_8_18/input/mol_1.xyz";
char *file2_name="/home/vanka/siddharth/ruchi/code_file_v1.0.5_with_counter_ions_8_8_18/input/mol_1.xyz";
int file_limit=142856; //file break point
int ring_atom=72;  //number of ring atoms
int track_start_atom_no=72;
int track_end_atom_no=207;
char *charge_file_name="/home/vanka/siddharth/ruchi/dm_test/charge_dm_test";


//directionality parameters
int start_frame=0;
int end_frame=142856;
int step_size=1;       //********The old "inc" variable
double epsi=37.6;
double K=23070.708216/37.6;     //in pN
//energy
double amu=1.6605e-27;
double angstorm=1e-10;
double femto=1e-15;
double timeStep=0.5;   //in femtoSeconds


//****notes*****//
// 1) Changed #atoms(atm_v) in shiftOrigin
// 2) Frame 2 should be used
// 3) variable "inc" is now called "step_size"
