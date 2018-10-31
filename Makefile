CC      =	gcc
CFLAGS  =	-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall  
MINIMAP2_DIR = ./minimap2
MINMAP2LIB  =   $(MINIMAP2_DIR)/libminimap2.a
PYLIB   =   -lpython2.7
LIB     =	$(MINMAP2LIB) -lm -lz -lpthread $(PYLIB)
PY_DIR = /usr/include/python2.7
INCLUDE = -I $(MINIMAP2_DIR) -I $(PY_DIR)
#TODO different -I for different .c

BIN_DIR =	./bin
SRC_DIR =   ./src

SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/miniTandem

GDB_DEBUG   =   $(BIN_DIR)/gdb_miniTandem
DMARCRO 	=	-D __DEBUG__

# dependencies
.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all:		$(MINMAP2LIB) $(BIN) 
miniTandem: $(BIN)
gdb_miniTandem: $(SOURCE) $(GDB_DEBUG) 



$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(OBJS) -o $@ $(LIB)


$(GDB_DEBUG):
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) $(INCLUDE) -o $@ $(LIB)

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
	rm -f $(SRC_DIR)/*.o $(GDB_DEBUG)
