CFLAGS = -O3 -fomit-frame-pointer -fkeep-inline-functions -finline-limit=50000 -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 -c
CLIBS = 
GCC = ~/MPI/bin/mpicc 
all: ExactMP

#run: ~/mpich/bin/mpirun $(NUMPROC) $(PROC) $(BINFILE) $(DATAFILE)

ExactMP: main.o preProc.o uTaskGen.o uCompCost.o uDataStr.o addOrder.o frameWork.o uBnB.o aggGreedy.o uPLB.o tbr.o eckDay.o nj.o randOpt.o oBnB.o oDataStr.o oPLB.o oTaskGen.o oDispTree.o opBnB.o
	$(GCC) main.o preProc.o uTaskGen.o uCompCost.o uDataStr.o addOrder.o frameWork.o uBnB.o aggGreedy.o uPLB.o tbr.o eckDay.o randOpt.o nj.o oBnB.o oDataStr.o oPLB.o oTaskGen.o opBnB.o oDispTree.o -lm -o ExactMP
	rm -rf *.o

main.o: main.c
	$(GCC) $(CFLAGS) main.c $(CLIBS) -o main.o

preProc.o: preProc.c
	$(GCC) $(CFLAGS) preProc.c $(CLIBS) -o preProc.o

uTaskGen.o: uTaskGen.c
	$(GCC) $(CFLAGS) uTaskGen.c $(CLIBS) -o uTaskGen.o

uCompCost.o: uCompCost.c
	$(GCC) $(CFLAGS) uCompCost.c $(CLIBS) -o uCompCost.o

uDataStr.o: uDataStr.c
	$(GCC) $(CFLAGS) uDataStr.c $(CLIBS) -o uDataStr.o

addOrder.o: addOrder.c 
	$(GCC) $(CFLAGS) addOrder.c $(CLIBS) -o addOrder.o

frameWork.o: frameWork.c 
	$(GCC) $(CFLAGS) frameWork.c $(CLIBS) -o frameWork.o

uBnB.o: uBnB.c 
	$(GCC) $(CFLAGS) uBnB.c $(CLIBS) -o uBnB.o

aggGreedy.o: aggGreedy.c 
	$(GCC) $(CFLAGS) aggGreedy.c $(CLIBS) -o aggGreedy.o

uPLB.o: uPLB.c
	$(GCC) $(CFLAGS) uPLB.c $(CLIBS) -o uPLB.o

tbr.o: tbr.c
	$(GCC) $(CFLAGS) tbr.c $(CLIBS) -o tbr.o

eckDay.o: eckDay.c
	$(GCC) $(CFLAGS) eckDay.c $(CLIBS) -o eckDay.o

nj.o: nj.c
	$(GCC) $(CFLAGS) nj.c $(CLIBS) -o nj.o

randOpt.o: randOpt.c
	$(GCC) $(CFLAGS) randOpt.c $(CLIBS) -o randOpt.o

oTaskGen.o: oTaskGen.c
	$(GCC) $(CFLAGS) oTaskGen.c $(CLIBS) -o oTaskGen.o

oBnB.o: oBnB.c
	$(GCC) $(CFLAGS) oBnB.c $(CLIBS) -o oBnB.o

oPLB.o: oPLB.c
	$(GCC) $(CFLAGS) oPLB.c $(CLIBS) -o oPLB.o

oDataStr.o: oDataStr.c
	$(GCC) $(CFLAGS) oDataStr.c $(CLIBS) -o oDataStr.o

oDispTree.o: oDispTree.c
	$(GCC) $(CFLAGS) oDispTree.c $(CLIBS) -o oDispTree.o

opBnB.o: opBnB.c
	$(GCC) $(CFLAGS) opBnB.c $(CLIBS) -o opBnB.o

clean:
	rm -rf *.o ExactMP
