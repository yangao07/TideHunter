CC      =	gcc
CPP		=   g++
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
CPPFLAGS=   -std=c++11
DFLAGS  =	-g -Wall #-pg
#MINIMAP2_DIR = ./minimap2
#MINMAP2LIB  =   $(MINIMAP2_DIR)/libminimap2.a
#PYLIB   =   -lpython2.7
INCLUDE =
LIB     =	-lm -lz -lpthread
LIB_DIR =   
PY_DIR = /usr/include/python2.7

EDLIB_DIR = ./edlib
EDLIB_INCLUDE = ./edlib/include
EDLIB     = $(EDLIB_DIR)/edlib.o

SPOA_DIR  = ./spoa
SPOA_INCLUDE = ./spoa/include/spoa
SPOALIB_DIR = $(SPOA_DIR)/build/lib
SPOALIB   = $(SPOA_DIR)/build/lib/libspoa.a
LIB_DIR   += -L $(SPOALIB_DIR)
LIB       += -lspoa

#INCLUDE = -I $(MINIMAP2_DIR) -I $(PY_DIR)
#TODO different -I for different .c

BIN_DIR =	./bin
SRC_DIR =   ./src

CSOURCE   = $(wildcard $(SRC_DIR)/*.c)
CPPSOURCE = $(wildcard $(SRC_DIR)/*.cpp)

SOURCE    = $(CSOURCE) $(CPPSOURCE)
DSOURCE   = $(SOURCE) $(EDLIB_DIR)/src/edlib.cpp

OBJS    =  $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
OBJS    += $(EDLIB)

BIN     =	$(BIN_DIR)/miniTandem

GDB_DEBUG   =   $(BIN_DIR)/gdb_miniTandem
DMARCRO 	=	-D __DEBUG__

# dependencies
.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

.cpp.o:
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) $< -o $@


all:		$(MINMAP2LIB) $(BIN) 
miniTandem: $(BIN)
gdb_miniTandem: $(SOURCE) $(GDB_DEBUG) 


$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CPP) $(OBJS) -L $(SPOALIB_DIR) $(LIB) -o $@

$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h
	$(CPP) -c $< -I $(EDLIB_INCLUDE) -o $@

$(SRC_DIR)/edlib_align.o: $(SRC_DIR)/edlib_align.c $(SRC_DIR)/edlib_align.h 
	$(CC) -c $(CFLAGS) $< -I $(EDLIB_INCLUDE) -o $@

$(SPOALIB): 
	wd=$(PWD); \
	cd $(SPOA_DIR); mkdir build; cd build; \
	cmake -DCMAKE_BUILD_TYPE=Release ..;	\
	make; cd $(wd);

$(SRC_DIR)/spoa_align.o: $(SRC_DIR)/spoa_align.cpp $(SRC_DIR)/spoa_align.h $(SPOALIB)
	$(CPP) -c $< $(CPPFLAGS) -I $(SPOA_INCLUDE) -o $@


$(GDB_DEBUG): $(DSOURCE)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CPP) $(CPPFLAGS) $(DFLAGS) $(DSOURCE) $(DMARCRO) $(INCLUDE) -I $(EDLIB_INCLUDE) -I $(SPOA_INCLUDE) $(LIB_DIR) $(LIB) -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
	rm -f $(GDB_DEBUG)
