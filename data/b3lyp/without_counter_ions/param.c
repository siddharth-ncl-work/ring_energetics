//** file for all global variables **//

int ax=0;
char *file1_name="/home/vanka/siddharth/ruchi/input_data/b3lyp/without_counter_ions/mol.xyz";
char *file2_name="/home/vanka/siddharth/ruchi/input_data/b3lyp/without_counter_ions/mol.xyz";
int file_limit=20127; //file break point
int ring_atom=72;  //number of ring atoms
int track_start_atom_no=72;
int track_end_atom_no=207;
char *charge_file_name="/home/vanka/siddharth/ruchi/dm_test/charge_dm_test";
char *version="1.0.0";
char *outfile_name="output/energy_data_b3lyp_without_counter_ions_v1.0.0_9_10_18.csv";

//directionality parameters
int start_frame=0;
int end_frame=20127;
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
