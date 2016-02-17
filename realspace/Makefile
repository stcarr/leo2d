CC_OBJ = sheet.o hstruct.o sdata.o intralayer_coupling.o interlayer_coupling.o
MPI_OBJ = main.o locality.o

OBJS = $(CC_OBJ) $(MPI_OBJ)

CC_SRC = main.cpp sheet.cpp hstruct.cpp sdata.cpp intralayer_coupling.cpp
MPI_SRC = locality.cpp
SRC_DIR = src
OUT_DIR = obj
EXECUTABLE = main

CC = gcc
MPICC = mpicxx
#DEBUG = -g -O0
DEBUG = 

CFLAGS = -c
MPIFLAGS = -c
LFLAGS = -L/n/home12/stcarr/soft/petsc/linux-c-debug/lib/ -lpetsc 

all: $(addprefix $(OUT_DIR)/,$(OBJS)) $(EXECUTABLE)

$(EXECUTABLE) : $(addprefix $(OUT_DIR)/,$(OBJS))
	$(MPICC) $(DEBUG)  $(LFLAGS) $(addprefix $(OUT_DIR)/,$(OBJS)) -o $@ -L. -lm
 	
$(addprefix $(OUT_DIR)/,$(CC_OBJ)): COMP = $(CC)
$(addprefix $(OUT_DIR)/,$(MPI_OBJ)): COMP = $(MPICC)
$(addprefix $(OUT_DIR)/,$(OBJS)): TARGET = $@

$(addprefix $(OUT_DIR)/,$(OBJS)) : $(addprefix $(SRC_DIR)/,$(OBJS:.o=.cpp))
	$(COMP) -I/n/home12/stcarr/soft/petsc/linux-c-debug/include/  $(DEBUG) -c $(addprefix $(SRC_DIR)/,$(notdir $(TARGET:.o=.cpp))) -o $@

.PHONY: clean

clean:
	rm -f  $(EXECUTABLE) obj/*.o


