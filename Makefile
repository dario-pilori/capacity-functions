CC = gcc
CFLAGS = -O2 -march=native
LDFLAGS = -lm

all: qam_gmi pam_gmi

pam_gmi: pam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

qam_gmi: qam_gmi.o gausshermite_functions.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
