CC_OBJ = sheet.o hstruct.o sdata.o intralayer_coupling.o interlayer_coupling.o
MPI_OBJ = main.o locality.o

OBJS = $(CC_OBJ) $(MPI_OBJ)

CC_SRC = main.cpp sheet.cpp hstruct.cpp sdata.cpp intralayer_coupling.cpp
MPI_SRC = locality.cpp
SRC_DIR = src
OUT_DIR = obj
EXECUTABLE = main

# ODYSSEY FILT_DIR = /n/kaxiras_lab/stcarr/soft/FILTLAN
FILT_DIR = /n/kaxiras_lab/stcarr/soft/FILTLAN

CC = g++
MPICC = mpicxx
DEBUG = -g -O0
#DEBUG = 

CFLAGS = -c
MPIFLAGS = -c
#LFLAGS = -L/n/home12/stcarr/soft/FILTLAN/LIB -L/n/regal/kaxiras_lab/locality_code/soft/petsc/linux-c-debug/lib -L/n/regal/kaxiras_lab/locality_code/lib -lpetsc

#MKLFLAGS = -lmkl_lapack95_ilp64 -lblas -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
#MKLFLAGS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_intel_thread -liomp5
MKLFLAGS = -lmkl_rt -lpthread -liomp5

LFLAGS = -L$(FILT_DIR)/LIB -ldfiltlan -ldmatkit $(MKLFLAGS)
IFLAGS = -DUSE_NAMESPACE -I$(FILT_DIR)/INCLUDE -I$(FILT_DIR)/MATKIT/INCLUDE
#IFLAGS = -DUSE_NAMESPACE -I$(FILT_DIR)/INCLUDE -I$(FILT_DIR)/MATKIT/INCLUDE -I/n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/intel-mkl/11.0.0.079-fasrc01/include

all: $(addprefix $(OUT_DIR)/,$(OBJS)) $(EXECUTABLE)

$(EXECUTABLE) : $(addprefix $(OUT_DIR)/,$(OBJS))
	$(MPICC)  $(addprefix $(OUT_DIR)/,$(OBJS)) -o $@ $(LFLAGS)
	
$(addprefix $(OUT_DIR)/,$(CC_OBJ)): COMP = $(CC)
$(addprefix $(OUT_DIR)/,$(MPI_OBJ)): COMP = $(MPICC)
$(addprefix $(OUT_DIR)/,$(OBJS)): TARGET = $@

$(addprefix $(OUT_DIR)/,$(OBJS)) : $(addprefix $(SRC_DIR)/,$(OBJS:.o=.cpp))
	$(COMP) $(IFLAGS) $(DEBUG) -c $(addprefix $(SRC_DIR)/,$(notdir $(TARGET:.o=.cpp))) -o $@

.PHONY: clean

clean:
	rm -f  $(EXECUTABLE) obj/*.o


