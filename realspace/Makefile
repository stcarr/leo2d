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
DEBUG = -g -O0
#DEBUG = 

CFLAGS = -c
MPIFLAGS = -c
#LFLAGS = -L/n/home12/stcarr/soft/FILTLAN/LIB -L/n/regal/kaxiras_lab/locality_code/soft/petsc/linux-c-debug/lib -L/n/regal/kaxiras_lab/locality_code/lib -lpetsc
LFLAGS = -L/home/tfuser/software/FILTLAN/LIB -ldfiltlan -ldmatkit -lm -lblas -llapack
IFLAGS = -DUSE_NAMESPACE -I/home/tfuser/software/FILTLAN/INCLUDE -I/home/tfuser/software/FILTLAN/MATKIT/INCLUDE

all: $(addprefix $(OUT_DIR)/,$(OBJS)) $(EXECUTABLE)

$(EXECUTABLE) : $(addprefix $(OUT_DIR)/,$(OBJS))
	$(MPICC) $(DEBUG) $(addprefix $(OUT_DIR)/,$(OBJS)) -o $@ $(LFLAGS)
	
$(addprefix $(OUT_DIR)/,$(CC_OBJ)): COMP = $(CC)
$(addprefix $(OUT_DIR)/,$(MPI_OBJ)): COMP = $(MPICC)
$(addprefix $(OUT_DIR)/,$(OBJS)): TARGET = $@

$(addprefix $(OUT_DIR)/,$(OBJS)) : $(addprefix $(SRC_DIR)/,$(OBJS:.o=.cpp))
	$(COMP) $(IFLAGS) $(DEBUG) -c $(addprefix $(SRC_DIR)/,$(notdir $(TARGET:.o=.cpp))) -o $@

.PHONY: clean

clean:
	rm -f  $(EXECUTABLE) obj/*.o


