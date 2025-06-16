include ../Makefile.$(SSW_COMPILE_TYPE)
include ../Makefile.etc

CFLAGS		= -c -I$(IDL_DIR)/external
FFLAGS		= -c

# Shorthand 
.SUFFIXES: .$(SO_EXT)

comp_poly.o: comp_poly.c

comp_poly.$(SO_EXT): comp_poly.o

all: comp_poly.$(SO_EXT)

install: comp_poly.$(SO_EXT) comp_poly.dlm
	mkdir -p $(SSW_DLM_TOP)/$(SSW_BINARY_TYPE)
	cp -p $^ $(SSW_DLM_TOP)/$(SSW_BINARY_TYPE)/

.c.o:
	$(CC) $(CFLAGS) $(X_CFLAGS) $*.c

.o.$(SO_EXT):
	$(LD) $(X_LD_FLAGS) -o $@ $?
