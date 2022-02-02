# Makefile for tapwdosp
OPTION = -O3
FORT = ifort 
tapwdoss2: tapwdoss2.o tsp04.o $(HOME)/lib/gaussj.o
	$(FORT) $(OPTION) tapwdoss2.o tsp04.o $(HOME)/lib/gaussj.o \
	-o tapwdoss2
prmdos.f:	prmdos.f
	cp prmdos.f prmdos.f
tapwdoss2.o: tapwdoss2.f prmdos.f
	$(FORT) $(OPTION) -c tapwdoss2.f
tapwdoss2.f: tapwdoss2.f
	ln -s tapwdoss2.f tapwdoss2.f
tsp04.o: tsp04.f prmtsp.f
	$(FORT) $(OPTION) -c tsp04.f
clean:
	rm *.o tapwdoss2
