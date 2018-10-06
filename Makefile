compile : main_rotation.c source/rotation.c source/shared_methods.c source/header_rotation.h source/error_analysis.c param.c source/code_verification.c source/translation.c source/dm_rotation.c
	gcc -Wall -ggdb3 main_rotation.c source/rotation.c source/shared_methods.c source/torque.c param.c source/error_analysis.c source/code_verification.c source/dm_rotation.c -lm -o output/a.out
	output/a.out
