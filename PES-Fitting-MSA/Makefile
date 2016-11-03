FC = ifort -r8 -O
LIBS = -mkl=sequential

%.o: %.f90
	$(FC) -c $<

OBJFIT := basis.o fit.o

OBJGET := basis.o gradient.o pes_shell.o getpot.o

fit.x : $(OBJFIT) 
	$(FC) -o $@ $^ $(LIBS)

getpot.x : $(OBJGET)
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.mod *.x	
