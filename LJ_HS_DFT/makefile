# comments
OBJS   =fmsa.o\
 ifmsa.o\
 fmt.o\
 fmt_a.o\
 ini.o\
 main.o\
 res.o\
 thermo.o
 
solvation: $(OBJS)
	mpif90 -O4 -o run $(OBJS)
clean:
	-rm *.o
