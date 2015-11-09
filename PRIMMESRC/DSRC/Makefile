#------------------------------------------------------------------------------
# Compile and add to the target library the core, double precision functions 
#------------------------------------------------------------------------------
include ../../Make_flags

# Include the primme.h and common_numerical.h header files
#
INCLUDE = -I../COMMONSRC

.PHONY: clean lib libd

CSOURCE = convergence_d.c correction_d.c primme_d.c init_d.c \
          inner_solve_d.c main_iter_d.c factorize_d.c numerical_d.c ortho_d.c \
	  restart_d.c locking_d.c solve_H_d.c update_projection_d.c update_W_d.c

COBJS = convergence_d.o correction_d.o primme_d.o init_d.o \
        inner_solve_d.o main_iter_d.o factorize_d.o numerical_d.o ortho_d.o \
	restart_d.o locking_d.o solve_H_d.o update_projection_d.o update_W_d.o 

convergence_d.o: convergence_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c convergence_d.c

correction_d.o: correction_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c correction_d.c

primme_d.o: primme_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c primme_d.c

init_d.o: init_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c init_d.c

inner_solve_d.o: inner_solve_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c inner_solve_d.c

main_iter_d.o: main_iter_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c main_iter_d.c

factorize_d.o: factorize_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c factorize_d.c

numerical_d.o: numerical_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c numerical_d.c

ortho_d.o: ortho_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c ortho_d.c

restart_d.o: restart_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c restart_d.c

locking_d.o: locking_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c locking_d.c

solve_H_d.o: solve_H_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c solve_H_d.c

update_projection_d.o: update_projection_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c update_projection_d.c

update_W_d.o: update_W_d.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c update_W_d.c

#
# Archive double precision object files in the full library
#
lib: $(COBJS)
	ar rv ../../$(LIBRARY) $(COBJS)
	$(RANLIB) ../../$(LIBRARY)
#
# Archive double precision object files in double-only library
#
libd: $(COBJS)
	ar rv ../../$(DLIBRARY) $(COBJS)
	$(RANLIB) ../../$(DLIBRARY)

#
# Clean up
#
clean:
	rm -f *.o a_d.out core

include dependencies
