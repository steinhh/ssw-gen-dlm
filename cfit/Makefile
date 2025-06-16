include ../Makefile.$(SSW_COMPILE_TYPE)

CFLAGS		= -c -I$(IDL_DIR)/external
FFLAGS          = -c

# Shorthand
E = $(SO_EXT)
.SUFFIXES: .$(E)

all: comp_poly.$(E)

install: comp_poly.$(E) comp_poly.dlm
	mkdir -p $(SSW_DLM_TOP)/$(SSW_BINARY_TYPE)
	cp -p $^ $(SSW_DLM_TOP)/$(SSW_BINARY_TYPE)/

.c.o:
	$(CC) $(CFLAGS) $(X_CFLAGS) $*.c

.o.$(E): 
	$(LD) $(X_LD_FLAGS) -o $@ $?

comp_poly.o: comp_poly.c

include ../Makefile.etc
