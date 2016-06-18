#-----------------------------------------------------------------
# Builds lib that includes both complex and double precision,
# or alternatively, one of the two precisions.
#-----------------------------------------------------------------
include ../Make_flags

libs: libd libz 
.PHONY: lib libs libd libz clean

lib: 
	@(\
	rm -f ../$(LIBRARY);\
	echo " Making common .o files"; \
	make -C COMMONSRC lib;\
	make -C SVDS/COMMONSRC lib;\
	echo " Making complex precision .o files"; \
	make -C ZSRC lib;\
	make -C SVDS/ZSRC lib;\
	echo " Making double precision .o files"; \
	make -C DSRC lib;\
	make -C SVDS/DSRC lib;\
	$(RANLIB) ../$(LIBRARY);\
	)

libd:
	@(\
	rm -f ../$(DLIBRARY);\
	echo " Making common .o files"; \
	make -C COMMONSRC libd;\
	make -C SVDS/COMMONSRC libd;\
	echo " Making double precision .o files"; \
	make -C DSRC libd;\
	make -C SVDS/DSRC libd;\
	$(RANLIB) ../$(DLIBRARY);\
	)

libz:
	@(\
	rm -f ../$(ZLIBRARY);\
	echo " Making common .o files"; \
	make -C COMMONSRC libz;\
	make -C SVDS/COMMONSRC libz;\
	echo " Making complex precision .o files"; \
	make -C ZSRC libz;\
	make -C SVDS/ZSRC libz;\
	$(RANLIB) ../$(ZLIBRARY);\
	)
#-----------------------------------------------------------------

clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " From PRIMME library directories"; \
	echo "   "COMMONSRC;\
	make -C COMMONSRC clean;\
	make -C SVDS/COMMONSRC clean;\
	echo "   "DSRC;\
	make -C DSRC clean;\
	make -C SVDS/DSRC clean;\
	echo "   "ZSRC;\
	make -C ZSRC clean;\
	make -C SVDS/ZSRC clean;\
	echo "--------------------------------------------------"; \
	)

