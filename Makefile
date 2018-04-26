CC = gcc
MEX = mex
CFLAGS = -O2 -fopenmp -fPIC
LDFLAGS = -lm 

all: qam_gmi pam_gmi qam_gmi_sweep qam_gmi_mex qam_llr_mex

qam_gmi_mex: gausshermite_functions.o
	$(MEX) -v CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_llr_mex: gausshermite_functions.o
	$(MEX) -v CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

pam_gmi: pam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi: qam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi_sweep: qam_gmi_sweep.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o pam_gmi qam_gmi qam_gmi_sweep *.mexa64
