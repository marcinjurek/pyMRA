include ../Makefile.in

.f.o:
	$(FC) $(FFLAGS) -fPIC -shared -c $<

OBJS =assmb.o \
btree2.o \
betree.o \
bfinit.o \
blkfct.o \
blkfc2.o \
blkslv.o \
chordr.o \
chlsup.o \
create.o \
dscal.o \
epost2.o \
etordr.o \
etpost.o \
etree.o \
fcnthn.o \
fsup1.o \
fsup2.o \
fnsplt.o \
fntsiz.o \
genmmd.o \
getrhs.o \
gtimer.o \
igathr.o \
inpnv.o \
invinv.o \
ldindx.o \
lstats.o \
mmdint.o \
mmdelm.o \
mmdupd.o \
mmdnum.o \
mmpy.o \
mmpy1.o \
mmpy2.o \
mmpy4.o \
mmpy8.o \
mmpyi.o \
ordmmd.o \
ordnat.o \
pchol.o \
sfinit.o \
smxpy1.o \
smxpy2.o \
smxpy4.o \
smxpy8.o \
symfct.o \
symfc2.o \
flo2ho.o \
ilo2ho.o \
dreadhb.o \
dsymmv.o \
selinv.o\
dumpL.o\
readmatrix.o\
readmatrixheader.o\
blkdinv.o

lib: $(OBJS)
	$(AR) libsupldlt.a $(OBJS)
	$(RANLIB) libsupldlt.a

clean:
	rm -f *.o libsupldlt.a
