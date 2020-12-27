CFLAGS = -O3 -fomit-frame-pointer -fkeep-inline-functions -finline-limit=50000 -fforce-addr -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4 -g -c
all: vc

vc: main.o preProc.o uTaskGen.o uCompCost.o uDataStr.o addOrder.o frameWork.o uBnB.o aggGreedy.o uPLB.o tbr.o eckDay.o nj.o randOpt.o oBnB.o oDataStr.o oPLB.o oTaskGen.o oDispTree.o
	gcc main.o preProc.o uTaskGen.o uCompCost.o uDataStr.o addOrder.o frameWork.o uBnB.o aggGreedy.o uPLB.o tbr.o eckDay.o randOpt.o nj.o oBnB.o oDataStr.o oPLB.o oTaskGen.o oDispTree.o -lm -o vc

main.o: main.c
	gcc $(CFLAGS) main.c -o main.o

preProc.o: preProc.c
	gcc $(CFLAGS) preProc.c -o preProc.o

uTaskGen.o: uTaskGen.c
	gcc $(CFLAGS) uTaskGen.c -o uTaskGen.o

uCompCost.o: uCompCost.c
	gcc $(CFLAGS) uCompCost.c -o uCompCost.o

uDataStr.o: uDataStr.c
	gcc $(CFLAGS) uDataStr.c -o uDataStr.o

addOrder.o: addOrder.c 
	gcc $(CFLAGS) addOrder.c -o addOrder.o

frameWork.o: frameWork.c 
	gcc $(CFLAGS) frameWork.c -o frameWork.o

uBnB.o: uBnB.c 
	gcc $(CFLAGS) uBnB.c -o uBnB.o

aggGreedy.o: aggGreedy.c 
	gcc $(CFLAGS) aggGreedy.c -o aggGreedy.o

uPLB.o: uPLB.c
	gcc $(CFLAGS) uPLB.c -o uPLB.o

tbr.o: tbr.c
	gcc $(CFLAGS) tbr.c -o tbr.o

eckDay.o: eckDay.c
	gcc $(CFLAGS) eckDay.c -o eckDay.o

nj.o: nj.c
	gcc $(CFLAGS) nj.c -o nj.o

randOpt.o: randOpt.c
	gcc $(CFLAGS) randOpt.c -o randOpt.o

oTaskGen.o: oTaskGen.c
	gcc $(CFLAGS) oTaskGen.c -o oTaskGen.o

oBnB.o: oBnB.c
	gcc $(CFLAGS) oBnB.c -o oBnB.o

oPLB.o: oPLB.c
	gcc $(CFLAGS) oPLB.c -o oPLB.o

oDataStr.o: oDataStr.c
	gcc $(CFLAGS) oDataStr.c -o oDataStr.o

oDispTree.o: oDispTree.c
	gcc $(CFLAGS) oDispTree.c -o oDispTree.o

clean:
	rm -f *.o vc
