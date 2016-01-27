OBJS = main.o sheet.o heterostructure.o
CC = gcc
DEBUG = -g
CFLAGS = -c
LFLAGS = 

main : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o main -lm -lstdc++

main.o : main.cpp sheet.hpp heterostructure.hpp
	$(CC) $(CFLAGS) main.cpp

sheet.o : sheet.hpp sheet.cpp
	$(CC) $(CFLAGS) sheet.cpp

heterostructure.o: heterostructure.hpp heterostructure.cpp sheet.hpp
	$(CC) $(CFLAGS) heterostructure.cpp

clean:
	rm -f *.o main


