compile:
gcc Source.c mmio.h mmio.compile

run:
./a.out 220.mtx

Change the state of these macros if required:
To enable filling out graph and running lexp/lexm on it to see if we get fillin edges
#define COMPARE_FILLINS 1
To allow lexm to run and generate an order and then feed it to GEM
#define ALLOWLEXM 1
To allow lexp to run and generate an order and then feed it to GEM
#define ALLOWLEXP 1

plotgraphs.m is a matlab script to plot the data.

scrrenshots and output files are included

Source_onlyGEM.c contains the GEM code. It does not compile. It just shows the relevant GEM code