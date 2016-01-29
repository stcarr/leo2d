OBJS = main.o sheet.o hstruct.o locality.o sdata.o
CC = gcc
MPICC = mpicxx
DEBUG = -g
CFLAGS = -c
MPIFLAGS = -c
LFLAGS = 

main : $(OBJS)
	$(MPICC) $(LFLAGS) $(OBJS) -o main -lm -lstdc++

main.o : main.cpp sheet.h hstruct.h sdata.h locality.h
	$(CC) $(CFLAGS) main.cpp

sheet.o : sheet.h sheet.cpp sdata.h
	$(CC) $(CFLAGS) sheet.cpp

hstruct.o: hstruct.h hstruct.cpp sheet.h sdata.h
	$(CC) $(CFLAGS) hstruct.cpp

locality.o: locality.h locality.cpp sheet.h hstruct.h sdata.h
	$(MPICC) $(MPIFLAGS) locality.cpp

sdata.o: sdata.h sdata.cpp
	$(CC) $(CFLAGS) sdata.cpp

clean:
	rm -f *.o main


