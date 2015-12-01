#------------------------------------------------------------------------------
# Compile and add to the target library the core, double precision functions 
#------------------------------------------------------------------------------
include ../../../Make_flags

# Include the primme_svds.h, primme.h and numerical_d.h header files
#
INCLUDE = -I../COMMONSRC -I../../COMMONSRC -I../../DSRC

.PHONY: clean lib libd

CSOURCE = primme_svds_d.c 

COBJS = primme_svds_d.o

primme_svds_d.o: primme_svds_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_svds_d.c


#
# Archive double precision object files in the full library
#
lib: $(COBJS)
	ar rv ../../../$(LIBRARY) $(COBJS)
	$(RANLIB) ../../../$(LIBRARY)
#
# Archive double precision object files in double-only library
#
libd: $(COBJS)
	ar rv ../../../$(DLIBRARY) $(COBJS)
	$(RANLIB) ../../../$(DLIBRARY)

#
# Clean up
#
clean:
	rm -f *.o a_d.out core

include dependencies
