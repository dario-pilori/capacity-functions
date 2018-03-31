CC = gcc
CFLAGS = -O2 -fopenmp
LDFLAGS = -lm 

all: qam_gmi pam_gmi qam_gmi_sweep

pam_gmi: pam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi: qam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi_sweep: qam_gmi_sweep.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
