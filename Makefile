CC      =   gcc
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-function -Wno-misleading-indentation

# for debug
ifneq ($(debug),)
	DFLAGS   = -D __DEBUG__
endif

# for gdb
ifneq ($(gdb),)
	CFLAGS   = -g -Wall ${DFLAGS}
else
	CFLAGS   += ${DFLAGS}
endif

# for gprof
ifneq ($(pg),)
	PG_FLAG  =   -pg
	CFLAGS  +=   -pg
endif

CPPFLAGS =  $(CFLAGS)

#TODO different -I for different .c
BIN_DIR =	./bin
SRC_DIR =   ./src
LIB_DIR =   ./lib
LIB           = -lm -lz -lpthread

INCLUDE       =

ifneq ($(BUILD_PREFIX),)
	INCLUDE_FLAG = "INCLUDE=-I$(BUILD_PREFIX)/include"
	LIB_FLAG = "-L$(BUILD_PREFIX)/lib"
endif


# edlib
EDLIB_DIR     = ./edlib
EDLIB_INCLUDE = ./edlib/include
EDLIB         = $(EDLIB_DIR)/src/edlib.o
#INCLUDE      += -I$(EDLIB_INCLUDE)

# SIMD
FLAG_SSE4       = -msse4
FLAG_AVX2       = -mavx2
FLAG_AVX512F    = -mavx512f
FLAG_AVX512BW   = -mavx512bw

# abPOA
ABPOA_DIR = ./abPOA
ABPOA_INCLUDE = $(ABPOA_DIR)/include
ABPOALIB      = ./lib/libabpoa.a
INCLUDE      += -I$(ABPOA_INCLUDE)
ABPOA_SIMD_FLAG =

# ksw2
# TODO flag with sse
KSW2_DIR     = ./ksw2
KSW2_INCLUDE = ./ksw2
#INCLUDE     += -I$(KSW2_INCLUDE)
KSW2_SIMD_FLAG = -msse4.1

ifneq ($(sse2),)
	KSW2_SIMD_FLAG = -msse2 -mno-sse4.1
	ABPOA_SIMD_FLAG = "sse2=1"
else ifneq ($(sse41),)
	KSW2_SIMD_FLAG = -msse4.1
	ABPOA_SIMD_FLAG = "sse41=1"
else ifneq ($(sse4),)
	KSW2_SIMD_FLAG = -msse4.1
	ABPOA_SIMD_FLAG = "sse41=1"
else ifneq ($(avx2),)
	KSW2_SIMD_FLAG = -msse4.1
	ABPOA_SIMD_FLAG = "avx2=1"
else ifneq ($(avx512f),)
	KSW2_SIMD_FLAG = -msse4.1
	ABPOA_SIMD_FLAG = "avx512f=1"
else ifneq ($(avx512bw),)
	KSW2_SIMD_FLAG = -msse4.1
	ABPOA_SIMD_FLAG = "avx512bw=1"
endif




CSOURCE    = $(wildcard $(SRC_DIR)/*.c)
CSOURCE   += $(wildcard $(KSW2_DIR)/*.c)
CPPSOURCE  = $(wildcard $(SRC_DIR)/*.cpp)
CPPSOURCE += $(EDLIB_DIR)/src/edlib.cpp

SOURCE  = $(CSOURCE) $(CPPSOURCE)
DSOURCE = $(SOURCE) 

OBJS    = $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
#OBJS   += $(EDLIB)

BIN     =	$(BIN_DIR)/TideHunter
ifneq ($(gdb),)
	BIN     =	$(BIN_DIR)/gdb_TideHunter
endif

# dependencies
.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

.cpp.o:
	$(CXX) -c $(CPPFLAGS) $(INCLUDE) $< -o $@

all:		$(MINMAP2LIB) $(BIN) 
TideHunter: $(BIN)

$(BIN): $(OBJS) $(ABPOALIB) Makefile
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CXX) $(CFLAGS) $(OBJS) $(ABPOALIB) $(INCLUDE) $(LIB) -o $@ $(PG_FLAG)

# edlib
$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h
	$(CXX) $(CFLAGS) -c $< $(INCLUDE) -I $(EDLIB_INCLUDE) -o $@

$(SRC_DIR)/edlib_align.o: $(SRC_DIR)/edlib_align.c $(SRC_DIR)/edlib_align.h 
	$(CC) -c $(CFLAGS) $< $(INCLUDE) -I $(EDLIB_INCLUDE) -o $@

# abPOA
$(ABPOALIB): 
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	cd $(ABPOA_DIR); \
	make libabpoa PREFIX=$(PWD) $(ABPOA_SIMD_FLAG) $(INCLUDE_FLAG) CFLAGS="$(LIB_FLAG) $(CFLAGS)"


# ksw2
$(KSW2_DIR)/ksw2_extz2_sse.o: $(KSW2_DIR)/ksw2_extz2_sse.c $(KSW2_DIR)/ksw2.h
	$(CC) -c $(CFLAGS) $(KSW2_SIMD_FLAG) $(INCLUDE) $< -o $@

$(KSW2_DIR)/ksw2_gg2_sse.o: $(KSW2_DIR)/ksw2_gg2_sse.c $(KSW2_DIR)/ksw2.h
	$(CC) -c $(CFLAGS) $(KSW2_SIMD_FLAG) $(INCLUDE) $< -o $@

#$(SRC_DIR)/abpoa_cons.o: $(SRC_DIR)/abpoa_cons.c $(ABPOA)
#	$(CC) -c $(CFLAGS) $< -I $(ABPOA_INCLUDE) $(INCLUDE) -o $@

# ksw2
$(SRC_DIR)/ksw2_align.o: $(SRC_DIR)/ksw2_align.c $(SRC_DIR)/ksw2_align.h
	$(CC) -c $(CFLAGS) $< -I $(KSW2_INCLUDE) $(INCLUDE) -o $@

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

clean_abPOA:
	rm -f $(ABPOALIB); \
	cd $(ABPOA_DIR); make clean; \

clean_ksw2:
	rm -f $(KSW2_DIR)/*.o

clean_edlib:
	rm -f $(EDLIB)

clean_all: clean clean_abPOA clean_ksw2 clean_edlib 

$(SRC_DIR)/abpoa_cons.o: $(SRC_DIR)/tidehunter.h $(SRC_DIR)/seq.h $(ABPOA_INCLUDE)/abpoa.h $(ABPOA_INCLUDE)/simd_instruction.h 
$(SRC_DIR)/edlib_align.o: $(SRC_DIR)/edlib_align.h $(SRC_DIR)/utils.h
$(SRC_DIR)/gen_cons.o: $(SRC_DIR)/utils.h $(SRC_DIR)/tidehunter.h $(SRC_DIR)/edlib_align.h $(SRC_DIR)/abpoa_cons.h $(SRC_DIR)/ksw2_align.h $(ABPOA_INCLUDE)/abpoa.h
$(SRC_DIR)/ksw2_align.o: $(SRC_DIR)/utils.h
$(SRC_DIR)/main.o: $(SRC_DIR)/utils.h $(SRC_DIR)/tidehunter.h $(SRC_DIR)/kseq.h $(SRC_DIR)/abpoa_cons.h $(ABPOA_INCLUDE)/abpoa.h
$(SRC_DIR)/partition.o: $(SRC_DIR)/utils.h $(SRC_DIR)/edlib_align.h $(SRC_DIR)/abpoa_cons.h $(SRC_DIR)/ksw2_align.h $(SRC_DIR)/tandem_chain.h
$(SRC_DIR)/seq.o: $(SRC_DIR)/utils.h $(SRC_DIR)/seq.h $(SRC_DIR)/kseq.h
$(SRC_DIR)/tandem_chain.o: $(SRC_DIR)/utils.h $(SRC_DIR)/tandem_chain.h $(SRC_DIR)/tandem_hit.h $(SRC_DIR)/tidehunter.h
$(SRC_DIR)/tandem_hit.o: $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h $(SRC_DIR)/tandem_hit.h $(SRC_DIR)/tidehunter.h
$(SRC_DIR)/tidehunter.o: $(SRC_DIR)/utils.h $(SRC_DIR)/tidehunter.h $(SRC_DIR)/tandem_chain.h $(SRC_DIR)/tandem_hit.h $(SRC_DIR)/partition.h $(SRC_DIR)/gen_cons.h $(SRC_DIR)/seq.h
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h
