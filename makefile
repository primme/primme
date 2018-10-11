#-----------------------------------------------------------------
# make 
#   lib	      makes the libprimme.a library
#   solib     makes the libprimme.so library
#   matlab    make libprimme.a compatible with MATLAB and the
#             module for MATLAB
#   octave    make libprimme.a and the Octave module
#   python    make libprimme.a and the python module
#   R_install install PRIMME's R interface
#   deps      update header dependencies
#   clean     removes all *.o files
#   clean_lib remove all files in lib
#   test      build and execute simple examples
#-----------------------------------------------------------------
include Make_flags

.PHONY: lib clean test all_tests check_style matlab octave python python_install R_install tags deps

#------------------------ Libraries ------------------------------
# Making the PRIMME library
# Includes float, double and complex counterparts
lib:
	@$(MAKE) -C src ../lib/$(LIBRARY)

solib:
	@$(MAKE) -C src ../lib/$(SONAMELIBRARY)
ifneq ($(SOLIBRARY),$(SONAMELIBRARY))
	@cd lib; ln -fs $(SONAMELIBRARY) $(SOLIBRARY)
endif
ifneq ($(MAJORVERSION),)
	@cd lib; ln -fs $(SONAMELIBRARY) $(SOLIBRARY).$(MAJORVERSION)
endif


clean: 
	@$(MAKE) -C src clean

clean_lib:
	@rm -f lib/*

test:
	@echo "------------------------------------------------";
	@echo " Test C examples                                ";
	@echo "------------------------------------------------";
	@$(MAKE) -C examples veryclean test_examples_C USE_PETSC=no

all_tests:
	@$(MAKE) -C examples veryclean test_examples;
	@$(MAKE) -C tests veryclean all_tests

matlab: clean clean_lib
	@$(MAKE) lib CFLAGS="${CFLAGS} -DPRIMME_BLASINT_SIZE=64"
	@$(MAKE) -C Matlab matlab

octave: clean clean_lib lib
	@$(MAKE) -C Matlab octave

python: clean clean_lib lib
	@$(MAKE) -C Python clean all

python_install: python
	@$(MAKE) -C Python install

R_install:
	@$(MAKE) -C R install

install: solib
	install -d $(includedir)
	install -m 644 include/primme.h include/primme_eigs.h \
		include/primme_eigs_f77.h include/primme_f77.h \
		include/primme_svds.h include/primme_svds_f77.h \
		$(includedir)
	install -d $(libdir)
	install -m 644 lib/$(SONAMELIBRARY) $(libdir)
ifneq ($(SOLIBRARY),$(SONAMELIBRARY))
	@cd $(libdir); ln -fs $(SONAMELIBRARY) $(SOLIBRARY)
endif
ifneq ($(MAJORVERSION),)
	@cd $(libdir); ln -fs $(SONAMELIBRARY) $(SOLIBRARY).$(MAJORVERSION)
endif

uninstall:
	rm -f $(libdir)/$(SONAMELIBRARY) $(libdir)/$(SOLIBRARY)
	rm -f $(includedir)/primme.h $(includedir)/primme_eigs.h \
		$(includedir)/primme_eigs_f77.h $(includedir)/primme_f77.h \
		$(includedir)/primme_svds.h $(includedir)/primme_svds_f77.h

deps:
	@$(MAKE) -C src deps
	@$(MAKE) -C src auto_headers

check_style:
	@( grep '	' -R . --include='*.[chfmF]' && echo "Please don't use tabs!" ) || true

tags:
	@$(MAKE) -C src ../tags

.NOTPARALLEL:
