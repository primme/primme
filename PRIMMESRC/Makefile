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
	echo " Making complex precision .o files"; \
	make -C ZSRC lib;\
	echo " Making double precision .o files"; \
	make -C DSRC lib;\
	$(RANLIB) ../$(LIBRARY);\
	)

libd:
	@(\
	rm -f ../$(DLIBRARY);\
	echo " Making common .o files"; \
	make -C COMMONSRC libd;\
	echo " Making double precision .o files"; \
	make -C DSRC libd;\
	$(RANLIB) ../$(DLIBRARY);\
	)

libz:
	@(\
	rm -f ../$(ZLIBRARY);\
	echo " Making common .o files"; \
	make -C COMMONSRC libz;\
	echo " Making complex precision .o files"; \
	make -C ZSRC libz;\
	$(RANLIB) ../$(ZLIBRARY);\
	)
#-----------------------------------------------------------------

clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " From PRIMME library directories"; \
	echo "   "COMMONSRC;\
	make -C COMMONSRC clean;\
	echo "   "DSRC;\
	make -C DSRC clean;\
	echo "   "ZSRC;\
	make -C ZSRC clean;\
	echo "--------------------------------------------------"; \
	)

