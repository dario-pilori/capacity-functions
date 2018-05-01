CC = gcc
MEX = mex
CFLAGS = -O2 -fopenmp -fPIC
LDFLAGS = -lm 

all: qam_gmi_mex qam_llr_mex calculate_pbit_mex qam_mi_montecarlo_mex

calculate_pbit_mex:
	$(MEX) CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_gmi_mex: capacity_functions.o
	$(MEX) CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_llr_mex: capacity_functions.o
	$(MEX) CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

pam_llr_mex: capacity_functions.o
	$(MEX) CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_mi_montecarlo_mex: capacity_functions.o
	$(MEX) CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

pam_gmi: pam_gmi.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi: qam_gmi.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi_sweep: qam_gmi_sweep.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o pam_gmi qam_gmi qam_gmi_sweep *.mexa64
