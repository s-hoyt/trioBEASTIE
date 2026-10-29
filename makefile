CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -I$(GSLDIR)/include -w
LDFLAGS		= $(OPTIMIZE) -L$(CONDA_PREFIX)/lib -Wl,-rpath-link=$(CONDA_PREFIX)/lib -Wl,-rpath=$(CONDA_PREFIX)/lib
BOOM		= BOOM
OBJ		= obj
LIBS		= -LBOOM -lBOOM -lgsl -lm -lgslcblas
#--------------------------------------------------------
$(OBJ)/phase-trio.o:\
		phase-trio.C
	$(CC) $(CFLAGS) -o $(OBJ)/phase-trio.o -c \
		phase-trio.C
#---------------------------------------------------------
phase-trio: \
		$(OBJ)/phase-trio.o
	$(CC) $(LDFLAGS) -o phase-trio \
		$(OBJ)/phase-trio.o \
		$(LIBS)
#---------------------------------------------------------