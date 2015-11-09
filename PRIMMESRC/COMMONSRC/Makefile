#------------------------------------------------------------------------------
# Compiles the common functions that are used by both complex and 
# double precisions codes. This includes all the interface to PRIMME
# except for the dprimme(), zprimme() functions.
#------------------------------------------------------------------------------
include ../../Make_flags

.PHONY: clean lib libd libz

CSOURCE =  Complexz.c common_numerical.c errors.c primme_f77.c \
	primme_interface.c wtime.c

COBJS = Complexz.o common_numerical.o  errors.o  primme_f77.o  \
	primme_interface.o wtime.o

# 
# Compilation
#
Complexz.o: Complexz.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c Complexz.c

common_numerical.o: common_numerical.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c common_numerical.c

errors.o: errors.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c errors.c

primme_f77.o: primme_f77.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_f77.c

primme_interface.o: primme_interface.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_interface.c

wtime.o: wtime.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c wtime.c

#
# Archive common object files in the target full library
#
lib: $(COBJS)
	ar rv ../../$(LIBRARY) $(COBJS)
	$(RANLIB) ../../$(LIBRARY)
#
# Archive common object files in the target Double library
#
libd: $(COBJS)
	ar rv ../../$(DLIBRARY) $(COBJS)
	$(RANLIB) ../../$(DLIBRARY)
#
# Archive common object files in the target Complex library
#
libz: $(COBJS)
	ar rv ../../$(ZLIBRARY) $(COBJS)
	$(RANLIB) ../../$(ZLIBRARY)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a.out core

include dependencies
