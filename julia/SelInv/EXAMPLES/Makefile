include ../Makefile.in
include makefile_allo

SUPLDLT = ../LIB/libsupldlt.a 

LIBS = $(METIS) $(SUPLDLT) $(LAPACK) $(BLAS) -lm 

OBJS = C2Finterface.o 

all: lap2d5pt.x selinv.x

.c.o:
	$(CC) $(CFLAGS) -fPIC -shared -c $<

.f.o:
	${FC} ${FFLAGS} -c $<

lap2d5pt.x: lap2d5pt.o $(OBJS) $(SUPLDLT) $(METIS)
	$(LD) $(LDFLAGS) lap2d5pt.o $(OBJS) $(LIBS) -o $@

selinv.x: selinv.o $(OBJS) $(SUPLDLT) $(METIS) 
	$(LD) $(LDFLAGS) selinv.o $(OBJS) $(LIBS) -o $@


$(SUPLDLT):
	(cd ../LIB; make lib)

julia:
	gcc -g -fPIC -shared -c C2Finterface1.c
	gcc -g -fPIC -shared -c selinvJulia.c
	gfortran  -g -DTIMING -fPIC -shared -llapack -lblas -lm C2Finterface1.o  selinvJulia.o $(allObjs) -o julia.so 

clean:
	rm -f *.o selinv.x lap2d5pt.x *.so
