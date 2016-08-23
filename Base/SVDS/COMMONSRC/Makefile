#

CDOTH = ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds.h  \
        ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.h \
        ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface_private.h \
        ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77_private.h \
        ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.h 

CSOURCE = ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.c \
          ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.c

all: $(CSOURCE) $(CDOTH) depends

#---------  PRIMMESRC/SVDS/COMMONSRC   .h files -------------------------

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds.h: primme_svds.h
	cp primme_svds.h ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds.h

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.h: primme_svds_interface.h
	cp primme_svds_interface.h ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.h

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface_private.h: primme_svds_interface_private.h
	cp primme_svds_interface_private.h ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface_private.h

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77_private.h: primme_svds_f77_private.h
	cp primme_svds_f77_private.h ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77_private.h

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.h: primme_svds_f77.h
	cp primme_svds_f77.h ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.h

#---------  PRIMMESRC/SVDS/COMMONSRC   .c files -------------------------

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.c: primme_svds_interface.c
	cp primme_svds_interface.c ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_interface.c

../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.c: primme_svds_f77.c
	cp primme_svds_f77.c ../../../PRIMMESRC/SVDS/COMMONSRC/primme_svds_f77.c

#
# Make dependencies in the two directories
#
#
depends:
	@(\
	cd ../../../PRIMMESRC/SVDS/COMMONSRC/;\
	gcc -I../../../PRIMMESRC/COMMONSRC/ -MM *.c > dependencies ; \
	)
#
#  clean	- remove all object files
#
clean:
	echo 'Cleaning up PRIMMESRC/SVDS/ZSRC and PRIMMESRC/SVDS/DSRC directories'
	rm -f $(CSOURCE_Z) $(FSOURCE_Z) $(CDOTH_Z)
	rm -f $(CSOURCE_D) $(CDOTH_D )
	rm -f a.out core
