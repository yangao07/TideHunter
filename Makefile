CC       =	 gcc
CPP		 =   g++
CFLAGS   =	 -Wall -g -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
CPPFLAGS =   -std=c++11

# for grpof
#CFLAGS   =	 -g -Wall -pg
#CPPFLAGS =   -std=c++11 -pg
#PG_FLAG  =   -pg

# for debug
#CFLAGS   =	 -g -Wall -D __DEBUG__ 

# use abPOA, otherwise use spoa
CFLAGS  +=   -D __ABPOA__

INCLUDE  =
LIB      =	 -lm -lz -lpthread

# edlib
EDLIB_DIR     = ./edlib
EDLIB_INCLUDE = ./edlib/include
EDLIB         = $(EDLIB_DIR)/src/edlib.o
INCLUDE      += -I$(EDLIB_INCLUDE)

# abPOA
ABPOA_DIR = ./abPOA
ABPOA_INCLUDE = $(ABPOA_DIR)/include
ABPOA_LIB_DIR = $(ABPOA_DIR)/lib
ABPOALIB      = $(ABPOA_DIR)/lib/libabpoa.a
ABPOALIB_FLAG = -L$(ABPOA_LIB_DIR) -labpoa
INCLUDE      += -I$(ABPOA_INCLUDE)

# ksw2
KSW2_DIR     = ./ksw2
KSW2_INCLUDE = ./ksw2
INCLUDE     += -I$(KSW2_INCLUDE)

# spoa
SPOA_DIR     = ./spoa
SPOA_INCLUDE = ./spoa/include/spoa
SPOALIB_DIR  = $(SPOA_DIR)/build/lib
SPOALIB      = $(SPOA_DIR)/build/lib/libspoa.a
SPOALIB_FLAG = -L$(SPOALIB_DIR) -lspoa
INCLUDE     += -I$(SPOA_INCLUDE)

#TODO different -I for different .c

BIN_DIR =	./bin
SRC_DIR =   ./src

CSOURCE    = $(wildcard $(SRC_DIR)/*.c)
CSOURCE   += $(wildcard $(KSW2_DIR)/*.c)
CPPSOURCE  = $(wildcard $(SRC_DIR)/*.cpp)
CPPSOURCE += $(EDLIB_DIR)/src/edlib.cpp

SOURCE  = $(CSOURCE) $(CPPSOURCE)
DSOURCE = $(SOURCE) 

OBJS    = $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
#OBJS   += $(EDLIB)

BIN     =	$(BIN_DIR)/miniTandem

#GDB_DEBUG   = $(BIN_DIR)/gdb_miniTandem
#ALL_GDB_DEBUG = $(BIN_DIR)/gdb_all_miniTandem
DMARCRO 	= -D __DEBUG__
ALL_HIT     = -D __ALL_HIT__ 

# dependencies
.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

.cpp.o:
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) $< -o $@


all:		$(MINMAP2LIB) $(BIN) 
miniTandem: $(BIN)
#gdb_miniTandem: $(SOURCE) $(GDB_DEBUG) 
#gdb_all_miniTandem: $(SOURCE) $(ALL_GDB_DEBUG) 


$(BIN): $(OBJS) $(ABPOALIB) $(SPOALIB) Makefile
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CPP) $(OBJS) $(ABPOALIB_FLAG) $(SPOALIB_FLAG) $(LIB) -o $@ $(PG_FLAG)

# edlib
$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h
	$(CPP) -c $< -I $(EDLIB_INCLUDE) -o $@

$(SRC_DIR)/edlib_align.o: $(SRC_DIR)/edlib_align.c $(SRC_DIR)/edlib_align.h 
	$(CPP) -c $(CPPFLAGS) $< -I $(EDLIB_INCLUDE) -o $@

# abPOA
$(ABPOALIB): 
	wd=$(PWD); \
	cd $(ABPOA_DIR); \
	make simd_check; make libabpoa; \
	cd $(wd);

# spoa
$(SPOALIB): 
	wd=$(PWD); \
	cd $(SPOA_DIR); mkdir build; cd build; \
	cmake -DCMAKE_BUILD_TYPE=Release ..;	\
	make; \
	cd $(wd);

$(SRC_DIR)/abpoa_align.o: $(SRC_DIR)/abpoa_align.c $(ABPOA)
	$(CC) -c $(CFLAGS) $< -I $(ABPOA_INCLUDE) -o $@

# ksw2
$(SRC_DIR)/ksw2_align.o: $(SRC_DIR)/ksw2_align.c $(SRC_DIR)/ksw2_align.h
	$(CC) -c $(CFLAGS) $< -I $(KSW2_INCLUDE) -o $@

$(SRC_DIR)/spoa_align.o: $(SRC_DIR)/spoa_align.cpp $(SRC_DIR)/spoa_align.h $(SPOALIB)
	$(CPP) -c $< $(CPPFLAGS) -I $(SPOA_INCLUDE) -o $@


$(ALL_GDB_DEBUG): $(DSOURCE) $(SPOALIB) Makefile
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CPP) $(CPPFLAGS) $(DFLAGS) $(DSOURCE) $(DMARCRO) $(ALL_HIT) $(INCLUDE) $(ABPOALIB_FLAG) $(SPOALIB_FLAG) $(LIB) -o $@

$(GDB_DEBUG): $(DSOURCE) $(ABPOALIB) $(SPOALIB) Makefile
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CPP) $(CPPFLAGS) $(DFLAGS) $(DSOURCE) $(DMARCRO) $(INCLUDE) $(ABPOALIB_FLAG) $(SPOALIB_FLAG) $(LIB) -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
	rm -f $(GDB_DEBUG)



$(SRC_DIR)/seq.o: $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h
$(SRC_DIR)/self_chain.o: $(SRC_DIR)/mini_tandem.h $(SRC_DIR)/self_chain.h $(SRC_DIR)/edlib_align.h $(SRC_DIR)/spoa_align.h $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h 
$(SRC_DIR)/main.o: $(SRC_DIR)/utils.h $(SRC_DIR)/mini_tandem.h
$(SRC_DIR)/mini_tandem.o: $(SRC_DIR)/mini_tandem.h $(SRC_DIR)/self_chain.h $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h
$(SRC_DIR)/edlib_align.o: $(SRC_DIR)/edlib_align.h
$(SRC_DIR)/spoa_align.o: $(SRC_DIR)/utils.h $(SRC_DIR)/spoa_align.h
$(SRC_DIR)/ksw2_align.o: $(SRC_DIR)/utils.h
