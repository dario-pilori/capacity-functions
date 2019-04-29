CC = gcc
MEX = mex
CFLAGS = -O2 -fopenmp -march=native
MEX_CFLAGS = -O2 -fopenmp -fPIC -march=native
LDFLAGS = -lm -fopenmp

all: qam_gmi_mex qam_llr_mex qam_mi_montecarlo_mex pam_mi_montecarlo_mex qam_llr_pn_mex qam_llr_pn_maxlog_mex qam_symbllr_mex qam_llr_maxlog_mex

qam_llr_pn_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_llr_pn_maxlog_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

calculate_pbit_mex:
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_gmi_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_llr_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_llr_maxlog_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

pam_llr_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_mi_montecarlo_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

pam_mi_montecarlo_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

qam_symbllr_mex: capacity_functions.o
	$(MEX) CFLAGS="$(MEX_CFLAGS)" LDFLAGS="$(LDFLAGS)" -R2018a $< $@.c

example_pam_gmi: example_pam_gmi.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

example_qam_gmi: example_qam_gmi.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

example_qam_gmi_sweep: example_qam_gmi_sweep.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

example_qam_gmi_lambda_sweep: example_qam_gmi_lambda_sweep.o capacity_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o pam_gmi qam_gmi qam_gmi_sweep *.mexa64
