CC = gcc
CFLAGS = -O
FC = gfortran
FFLAGS = -O
LIBS = -lm
PROGS=admix count hetfst hwe hwtest inform ldtest

all: $(PROGS)

admix.o: admix.c
	$(CC) $(CFLAGS) -c admix.c

inform.o: inform.c
	$(CC) $(CFLAGS) -c inform.c

dgedifa.o: dgedifa.f
	$(FC) $(FFLAGS) -c dgedifa.f

admix: admix.o dgedifa.o
	$(FC) -o admix admix.o dgedifa.o $(LIBS)

inform: inform.o dgedifa.o
	$(FC) -o inform inform.o dgedifa.o $(LIBS)

count: count.c
	$(CC) $(CFLAGS) -o count count.c $(LIBS)

hetfst: hetfst.c
	$(CC) $(CFLAGS) -o hetfst hetfst.c $(LIBS)

hwtest.o: hwtest.c
	$(CC) $(CFLAGS) -c hwtest.c

ldtest.o: ldtest.c
	$(CC) $(CFLAGS) -c ldtest.c

pchis.o: pchis.f
	$(FC) $(FFLAGS) -c pchis.f

hwtest: hwtest.o pchis.o
	$(FC) -o hwtest hwtest.o pchis.o $(LIBS)

ldtest: ldtest.o pchis.o
	$(FC) -o ldtest ldtest.o pchis.o $(LIBS)

hwe:
	cd HWE-src ; make hwe ; make clean ; mv hwe ..

clean:
	rm -f *.o

distclean: clean
	rm -f $(PROGS)
