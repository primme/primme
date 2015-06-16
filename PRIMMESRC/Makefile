#-----------------------------------------------------------------
# Builds lib that includes both complex and double precision,
# or alternatively, one of the two precisions.
#-----------------------------------------------------------------
include ../Make_flags

libs: libd libz 
.PHONY: lib libs libd libz clean

lib: 
	@(\
	rm -f $(LIBRARY);\
	echo " Making common .o files"; \
	cd $(COMMONDIR) ;\
	make lib;\
	echo " Making complex precision .o files"; \
	cd $(ZPRIMMEDIR) ;\
	make lib;\
	echo " Making double precision .o files"; \
	cd $(DPRIMMEDIR) ;\
	make lib;\
	$(RANLIB) $(LIBRARY);\
	)

libd:
	@(\
	rm -f $(DLIBRARY);\
	echo " Making common .o files"; \
	cd $(COMMONDIR) ;\
	make libd;\
	echo " Making double precision .o files"; \
	cd $(DPRIMMEDIR) ;\
	make libd;\
	$(RANLIB) $(DLIBRARY);\
	)

libz:
	@(\
	rm -f $(ZLIBRARY);\
	echo " Making common .o files"; \
	cd $(COMMONDIR) ;\
	make libz;\
	echo " Making double precision .o files"; \
	cd $(ZPRIMMEDIR) ;\
	make libz;\
	$(RANLIB) $(ZLIBRARY);\
	)
#-----------------------------------------------------------------

clean: 
	@(\
	echo "--------------------------------------------------"; \
	echo " Cleaning .o, and archives from all directories   "; \
	echo "--------------------------------------------------"; \
	cd $(COMMONDIR) ;\
	make -f Makefile clean;\
	cd $(DPRIMMEDIR) ;\
	make -f Makefile clean;\
	cd $(ZPRIMMEDIR) ;\
	make -f Makefile clean;\
	)

