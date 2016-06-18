#------------------------------------------------------------------------------
# Compile and add to the target library the core, complex functions 
#------------------------------------------------------------------------------
include ../../../Make_flags

# Include the primme_svds.h, primme.h and numerical_Z.h header files
#
INCLUDE =  -I../COMMONSRC/ -I../../COMMONSRC/ -I../../ZSRC/

.PHONY: clean lib libz

CSOURCE = primme_svds_z.c 

COBJS = primme_svds_z.o

primme_svds_z.o: primme_svds_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_svds_z.c


#
# Archive complex precision object files in full library
#
lib: $(COBJS) $(FOBJS)
	ar rv ../../../$(LIBRARY) $(COBJS) $(FOBJS)
	$(RANLIB) ../../../$(LIBRARY)
#
# Archive complex precision object files in Complex-only library
#
libz: $(COBJS) $(FOBJS)
	ar rv ../../../$(ZLIBRARY) $(COBJS) $(FOBJS)
	$(RANLIB) ../../../$(ZLIBRARY)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a_z.out core

include dependencies
