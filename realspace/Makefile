CC_OBJ = sheet.o hstruct.o sdata.o intralayer_coupling.o interlayer_coupling.o spmatrix.o
MPI_OBJ = main.o locality.o

OBJS = $(CC_OBJ) $(MPI_OBJ)

CC_SRC = main.cpp sheet.cpp hstruct.cpp sdata.cpp intralayer_coupling.cpp spmatrix.cpp
MPI_SRC = locality.cpp
SRC_DIR = src
OUT_DIR = obj
EXECUTABLE = LEO2D




CC = g++
MPICC = mpicxx

DEBUG = -g -O0
#DEBUG = 

IFLAGS = -DUSE_MKL

MKLFLAGS = -lmkl_rt -lpthread -liomp5
#MKLFLAGS = 

LFLAGS = $(MKLFLAGS)

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


