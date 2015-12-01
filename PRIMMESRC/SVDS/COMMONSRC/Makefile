#------------------------------------------------------------------------------
# Compiles the common functions that are used by both complex and 
# double precisions codes. This includes all the interface to PRIMME
# except for the dprimme_svds(), zprimme_svds() functions.
#------------------------------------------------------------------------------
include ../../../Make_flags

# Include the primme.h  header files
INCLUDE = -I../../COMMONSRC

.PHONY: clean lib libd libz

CSOURCE =  primme_svds_interface.c

COBJS = primme_svds_interface.o 

# 
# Compilation
#

primme_svds_interface.o: primme_svds_interface.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_svds_interface.c

#
# Archive common object files in the target full library
#
lib: $(COBJS)
	ar rv ../../../$(LIBRARY) $(COBJS)
	$(RANLIB) ../../../$(LIBRARY)
#
# Archive common object files in the target Double library
#
libd: $(COBJS)
	ar rv ../../../$(DLIBRARY) $(COBJS)
	$(RANLIB) ../../../$(DLIBRARY)
#
# Archive common object files in the target Complex library
#
libz: $(COBJS)
	ar rv ../../../$(ZLIBRARY) $(COBJS)
	$(RANLIB) ../../../$(ZLIBRARY)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a.out core

include dependencies
