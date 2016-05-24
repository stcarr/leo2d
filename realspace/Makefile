CC_OBJ = sheet.o hstruct.o sdata.o intralayer_coupling.o interlayer_coupling.o spmatrix.o
MPI_OBJ = main.o locality.o

OBJS = $(CC_OBJ) $(MPI_OBJ)

CC_SRC = main.cpp sheet.cpp hstruct.cpp sdata.cpp intralayer_coupling.cpp spmatrix.cpp
MPI_SRC = locality.cpp
SRC_DIR = src
OUT_DIR = obj
EXECUTABLE = LEO2D

DEBUG = -g -O0

PLATFORM = ODYSSEY
#WFLAG = -DLINE_WORK
WFLAG = -DSQ_WORK


# Block for Odyssey
ifeq ($(PLATFORM),ODYSSEY)
CC = g++
MPICC = mpicxx
PFLAGS = -DUSE_MKL
LOCALFLAGS = -lmkl_rt -lpthread -liomp5
endif

# Block for VESTA/MIRA
ifeq ($(PLATFORM),ALCF)
CC = mpixlc_r
MPICC = mpixlcxx_r
PFLAGS = 
LOCALFLAGS =
#LOCALFLAGS = -lmkl_rt -lpthread -liomp5 -L/soft/libraries/alcf/current/xl/LAPACK/lib -L/soft/libraries/alcf/current/xl/BLAS/lib -llapack -lblas
endif

# Generic Block
ifeq ($(PLATFORM),DEFAULT)
CC = g++
MPICC = mpicxx
PFLAGS = 
LOCALFLAGS = 
endif


LFLAGS = $(LOCALFLAGS)

all: $(addprefix $(OUT_DIR)/,$(OBJS)) $(EXECUTABLE)

$(EXECUTABLE) : $(addprefix $(OUT_DIR)/,$(OBJS))
	$(MPICC)  $(addprefix $(OUT_DIR)/,$(OBJS)) -o $@ $(LFLAGS)
	
$(addprefix $(OUT_DIR)/,$(CC_OBJ)): COMP = $(CC)
$(addprefix $(OUT_DIR)/,$(MPI_OBJ)): COMP = $(MPICC)
$(addprefix $(OUT_DIR)/,$(OBJS)): TARGET = $@

$(addprefix $(OUT_DIR)/,$(OBJS)) : $(addprefix $(SRC_DIR)/,$(OBJS:.o=.cpp))
	$(COMP) $(PFLAGS) $(WFLAG) $(DEBUG) -c $(addprefix $(SRC_DIR)/,$(notdir $(TARGET:.o=.cpp))) -o $@

.PHONY: clean

clean:
	rm -f  $(EXECUTABLE) obj/*.o


