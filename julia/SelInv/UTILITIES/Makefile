include ../Makefile.in

convertmat.x: convertmat.o mmio.o
	$(FC) $(FFLAGS) $^ -o $@

convertmat.o: convertmat.f90
	$(FC) $(FFLAGS) $< -c

mmio.o: mmio.f
	$(FC) $(FFLAGS) $< -c

clean: 
	rm -rf *.o *.x
