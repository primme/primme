
include ../Make_flags

HEADERS := \
   include/common.h \
   include/template.h \
   include/template_types.h \
   include/template_undef.h \
   include/numerical.h \
   include/primme_interface.h \
   include/memman.h \
   eigs/common_eigs.h \
   eigs/template_normal.h \
   svds/primme_svds_interface.h

SOURCES := \
   eigs/auxiliary_eigs.c \
   eigs/auxiliary_eigs_normal.c \
   eigs/convergence.c \
   eigs/correction.c \
   eigs/factorize.c \
   eigs/init.c \
   eigs/inner_solve.c \
   eigs/main_iter.c \
   eigs/ortho.c \
   eigs/primme_c.c \
   eigs/primme_f77.c \
   eigs/primme_interface.c \
   eigs/restart.c \
   eigs/solve_projection.c \
   eigs/update_projection.c \
   eigs/update_W.c \
   linalg/auxiliary.c \
   linalg/blaslapack.c \
   linalg/magma_wrapper.c \
   linalg/memman.c \
   linalg/wtime.c \
   svds/primme_svds_c.c \
   svds/primme_svds_f77.c \
   svds/primme_svds_interface.c

OBJS := $(patsubst %.c,%.o,$(SOURCES))

INCLUDE += -I../include -Iinclude
INCLUDES := $(sort $(foreach dir,../include,$(wildcard $(dir)/*.h))) $(HEADERS)

#
# Generation of automatic headers
#

AUTOMATED_HEADERS_LINALG := \
   include/blaslapack.h \
   include/magma_wrapper.h \
   include/auxiliary.h
AUTOMATED_HEADERS_OTHERS := \
   eigs/auxiliary_eigs.h \
   eigs/auxiliary_eigs_normal.h \
   eigs/solve_projection.h \
   eigs/convergence.h \
   eigs/inner_solve.h \
   eigs/main_iter.h \
   eigs/ortho.h \
   eigs/primme_c.h \
   eigs/factorize.h \
   eigs/restart.h \
   eigs/update_W.h \
   eigs/correction.h \
   eigs/update_projection.h \
   eigs/init.h \
   svds/primme_svds_c.h

CPPFLAGS_ONLY_PREPROCESS ?= -E

$(AUTOMATED_HEADERS_LINALG): include/%.h : linalg/%.c tools/AUTO_HEADER
	cat tools/AUTO_HEADER > $@; \
	echo "#ifndef $(*F)_H" >> $@; \
	echo "#define $(*F)_H" >> $@; \
	$(CC) $(CPPFLAGS_ONLY_PREPROCESS) $(CFLAGS) $(DEFINES) $(INCLUDE) -DCHECK_TEMPLATE $< | $(PYTHON) tools/ctemplate >> $@; \
	echo "#endif" >> $@

$(AUTOMATED_HEADERS_OTHERS): %.h : %.c tools/AUTO_HEADER
	cat tools/AUTO_HEADER > $@; \
	echo "#ifndef $(*F)_H" >> $@; \
	echo "#define $(*F)_H" >> $@; \
	$(CC) $(CPPFLAGS_ONLY_PREPROCESS) $(CFLAGS) $(DEFINES) $(INCLUDE) -DCHECK_TEMPLATE $< | $(PYTHON) tools/ctemplate >> $@; \
	echo "#endif" >> $@

auto_headers: $(AUTOMATED_HEADERS_LINALG) $(AUTOMATED_HEADERS_OTHERS)

#
# Compilation
#

%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@

#
# Archive object files in the full library
#

UNAME ?= $(shell uname)

../lib/$(LIBRARY): $(OBJS)
	@mkdir -p ../lib
	@rm -f $@
ifeq ($(UNAME), Darwin)
	@libtool -static -o $@ $(OBJS) -no_warning_for_no_symbols
else
	@$(AR) r $@ $(OBJS)
	@$(RANLIB) $@
endif

ifeq ($(UNAME), Darwin)
../lib/$(SONAMELIBRARY): $(OBJS)
	@mkdir -p ../lib
	@$(CC) -shared -Wl,-install_name,libprimme.so.$(MAJORVERSION) -Wl,-undefined -Wl,dynamic_lookup $(OBJS) -o ../lib/$(SONAMELIBRARY) $(LDFLAGS)
else
../lib/$(SONAMELIBRARY): $(OBJS)
	@mkdir -p ../lib
	@$(CC) -shared -Wl,-soname,libprimme.so.$(MAJORVERSION) -o ../lib/$(SONAMELIBRARY) -Wl,--whole-archive $(OBJS) -Wl,--no-whole-archive $(LDFLAGS)
endif

#
# Generate CTAGS
#

../tags: $(SOURCES) $(HEADERS)
	@ctags --tag-relative=yes  -o $@.base $^;\
	cp $@.base $@.all;\
	for pre in R SH RH s c d z; do \
		sed "s/Sprimme/$${pre}primme/" $@.base >> $@.all;\
	done;\
	LC_COLLATE=C sort -d -u $@.all > $@;\
	rm -f $@.base $@.all

#
# Cleaning
#

clean:
	-@rm -f $(OBJS)

.DELETE_ON_ERROR: $(OBJS) 
.PHONY: auto_headers clean

#
# Dependencies
#

# NOTE: make may update deps and reload it automatically
deps: $(INCLUDES) $(SOURCES) $(AUTOMATED_HEADERS_LINALG) $(AUTOMATED_HEADERS_OTHERS)
	echo "# This file is generated automatically. Please don't modify" > deps
	
	# Added dependencies on headers for building the objects
	$(PYTHON) tools/dependencies $(INCLUDE) $(SOURCES) -e.o >> deps

include deps

