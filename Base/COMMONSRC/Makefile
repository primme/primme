#

CDOTH = ../../PRIMMESRC/COMMONSRC/common_numerical.h \
	../../PRIMMESRC/COMMONSRC/common_numerical_private.h \
	../../PRIMMESRC/COMMONSRC/Complexz.h \
	../../PRIMMESRC/COMMONSRC/const.h \
	../../PRIMMESRC/COMMONSRC/errors_private.h  \
	../../PRIMMESRC/COMMONSRC/primme.h \
	../../PRIMMESRC/COMMONSRC/primme_f77.h  \
	../../PRIMMESRC/COMMONSRC/primme_f77_private.h \
	../../PRIMMESRC/COMMONSRC/wtime.h

CSOURCE = ../../PRIMMESRC/COMMONSRC/common_numerical.c \
	../../PRIMMESRC/COMMONSRC/Complexz.c \
	../../PRIMMESRC/COMMONSRC/errors.c \
	../../PRIMMESRC/COMMONSRC/primme_f77.c \
	../../PRIMMESRC/COMMONSRC/primme_interface.c \
	../../PRIMMESRC/COMMONSRC/wtime.c

all: $(CSOURCE) $(CDOTH) depends

#---------  PRIMMESRC/COMMONSRC   .h files -------------------------

../../PRIMMESRC/COMMONSRC/common_numerical.h: common_numerical.h
	cp common_numerical.h ../../PRIMMESRC/COMMONSRC/common_numerical.h

../../PRIMMESRC/COMMONSRC/common_numerical_private.h: common_numerical_private.h
	cp common_numerical_private.h ../../PRIMMESRC/COMMONSRC/common_numerical_private.h

../../PRIMMESRC/COMMONSRC/Complexz.h: Complexz.h
	cp Complexz.h ../../PRIMMESRC/COMMONSRC/Complexz.h 

../../PRIMMESRC/COMMONSRC/const.h: const.h
	cp const.h ../../PRIMMESRC/COMMONSRC/const.h 

../../PRIMMESRC/COMMONSRC/errors_private.h: errors_private.h
	cp errors_private.h ../../PRIMMESRC/COMMONSRC/errors_private.h 

../../PRIMMESRC/COMMONSRC/primme.h: primme.h
	cp primme.h ../../PRIMMESRC/COMMONSRC/primme.h

../../PRIMMESRC/COMMONSRC/primme_f77.h: primme_f77.h
	cp primme_f77.h ../../PRIMMESRC/COMMONSRC/primme_f77.h

../../PRIMMESRC/COMMONSRC/primme_f77_private.h: primme_f77_private.h
	cp primme_f77_private.h ../../PRIMMESRC/COMMONSRC/primme_f77_private.h

../../PRIMMESRC/COMMONSRC/wtime.h: wtime.h
	cp wtime.h ../../PRIMMESRC/COMMONSRC/wtime.h

#---------  PRIMMESRC/COMMONSRC   .c files -------------------------

../../PRIMMESRC/COMMONSRC/common_numerical.c: common_numerical.c
	cp common_numerical.c ../../PRIMMESRC/COMMONSRC/common_numerical.c

../../PRIMMESRC/COMMONSRC/Complexz.c: Complexz.c
	cp Complexz.c ../../PRIMMESRC/COMMONSRC/Complexz.c

../../PRIMMESRC/COMMONSRC/errors.c: errors.c
	cp errors.c ../../PRIMMESRC/COMMONSRC/errors.c 

../../PRIMMESRC/COMMONSRC/primme_interface.c: primme_interface.c
	cp primme_interface.c ../../PRIMMESRC/COMMONSRC/primme_interface.c

../../PRIMMESRC/COMMONSRC/primme_f77.c: primme_f77.c
	cp primme_f77.c ../../PRIMMESRC/COMMONSRC/primme_f77.c

../../PRIMMESRC/COMMONSRC/wtime.c: wtime.c
	cp wtime.c ../../PRIMMESRC/COMMONSRC/wtime.c

#
# Make dependencies in the two directories
#
#
depends:
	@(\
	cd ../../PRIMMESRC/COMMONSRC/;\
	gcc -MM *.c > dependencies ; \
	)
#
#  clean	- remove all object files
#
clean:
	echo 'Cleaning up PRIMMESRC/ZSRC and PRIMMESRC/DSRC directories'
	rm -f $(CSOURCE_Z) $(FSOURCE_Z) $(CDOTH_Z)
	rm -f $(CSOURCE_D) $(CDOTH_D )
	rm -f a.out core
