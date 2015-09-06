#------------------------------------------------------------------------------
# Compile and add to the target library the core, complex functions 
#------------------------------------------------------------------------------
include ../../Make_flags

# Include the primme.h and common_numerical.h header files
#
INCLUDE =  -I../COMMONSRC/

.PHONY: clean lib libz

CSOURCE = convergence_z.c correction_z.c primme_z.c init_z.c \
          inner_solve_z.c main_iter_z.c factorize_z.c numerical_z.c ortho_z.c \
	  restart_z.c locking_z.c solve_H_z.c update_projection_z.c update_W_z.c

COBJS = convergence_z.o correction_z.o primme_z.o init_z.o \
        inner_solve_z.o main_iter_z.o factorize_z.o numerical_z.o ortho_z.o \
	restart_z.o locking_z.o solve_H_z.o update_projection_z.o update_W_z.o 

convergence_z.o: convergence_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c convergence_z.c

correction_z.o: correction_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c correction_z.c

primme_z.o: primme_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_z.c

init_z.o: init_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c init_z.c

inner_solve_z.o: inner_solve_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c inner_solve_z.c

main_iter_z.o: main_iter_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c main_iter_z.c

factorize_z.o: factorize_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c factorize_z.c

numerical_z.o: numerical_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c numerical_z.c

ortho_z.o: ortho_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c ortho_z.c

restart_z.o: restart_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c restart_z.c

locking_z.o: locking_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c locking_z.c

solve_H_z.o: solve_H_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c solve_H_z.c

update_projection_z.o: update_projection_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c update_projection_z.c

update_W_z.o: update_W_z.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c update_W_z.c

#
# Archive complex precision object files in full library
#
lib: $(COBJS) $(FOBJS)
	ar rv ../../$(LIBRARY) $(COBJS) $(FOBJS)
	$(RANLIB) ../../$(LIBRARY)
#
# Archive complex precision object files in Complex-only library
#
libz: $(COBJS) $(FOBJS)
	ar rv ../../$(ZLIBRARY) $(COBJS) $(FOBJS)
	$(RANLIB) ../../$(ZLIBRARY)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a_z.out core

include dependencies
