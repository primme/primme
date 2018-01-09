
include ../Make_flags

UNAME ?= $(shell uname)

DIRS := linalg eigs svds
SOURCES := $(sort $(foreach dir,$(DIRS),$(wildcard $(dir)/*.c)))

ARCHS := float floatcomplex double doublecomplex
OBJS := $(foreach arch,$(ARCHS),$(patsubst %.c,%$(arch).o,$(SOURCES)))

OTHER_HEADERS := \
   include/template.h \
   include/numerical.h \
   include/primme_interface.h

INCLUDE += -I../include -Iinclude
INCLUDES := $(sort $(foreach dir,../include,$(wildcard $(dir)/*.h))) $(OTHER_HEADERS)

#
# Generation of automatic headers
#

AUTOMATED_HEADERS_LINALG := \
   include/blaslapack.h \
   include/auxiliary.h
AUTOMATED_HEADERS_EIGS := \
   eigs/auxiliary_eigs.h \
   eigs/solve_projection.h \
   eigs/convergence.h \
   eigs/inner_solve.h \
   eigs/main_iter.h \
   eigs/ortho.h \
   eigs/globalsum.h \
   eigs/factorize.h \
   eigs/restart.h \
   eigs/update_W.h \
   eigs/correction.h \
   eigs/update_projection.h \
   eigs/init.h

CPPFLAGS_ONLY_PREPROCESS ?= -E

$(AUTOMATED_HEADERS_LINALG): include/%.h : linalg/%.c tools/AUTO_HEADER
	@cat tools/AUTO_HEADER > $@; \
	echo "#ifndef $(*F)_H" > $@.t; \
	echo "#define $(*F)_H" >> $@.t; \
	for arch in DOUBLE DOUBLECOMPLEX FLOAT FLOATCOMPLEX ; do \
	   $(CC) $(CPPFLAGS_ONLY_PREPROCESS) $(CFLAGS) $(DEFINES) $(INCLUDE) -DCHECK_TEMPLATE -DUSE_$$arch $< | tools/ctemplate >> $@.t; \
	done; \
	cat $@.t >> $@; rm -f $@.t; \
	echo "#endif" >> $@

$(AUTOMATED_HEADERS_EIGS): %.h : %.c tools/AUTO_HEADER
	@cat tools/AUTO_HEADER > $@; \
	echo "#ifndef $(*F)_H" > $@.t; \
	echo "#define $(*F)_H" >> $@.t; \
	for arch in DOUBLE DOUBLECOMPLEX FLOAT FLOATCOMPLEX ; do \
	   $(CC) $(CPPFLAGS_ONLY_PREPROCESS) $(CFLAGS) $(DEFINES) $(INCLUDE) -DCHECK_TEMPLATE -DUSE_$$arch $< | tools/ctemplate >> $@.t; \
	done; \
	cat $@.t >> $@; rm -f $@.t; \
	echo "#endif" >> $@

auto_headers: $(AUTOMATED_HEADERS_LINALG) $(AUTOMATED_HEADERS_EIGS)

#
# Compilation
#

%double.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_DOUBLE -c $< -o $@

%doublecomplex.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_DOUBLECOMPLEX -c $< -o $@

%float.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_FLOAT -c $< -o $@

%floatcomplex.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_FLOATCOMPLEX -c $< -o $@

#
# Archive object files in the full library
#

../lib/$(LIBRARY): $(OBJS)
	@mkdir -p ../lib
ifeq ($(UNAME), Darwin)
	@libtool -static -o $@ $(OBJS) -no_warning_for_no_symbols
else
	@$(AR) r $@ $(OBJS)
	@$(RANLIB) $@
endif

ifeq ($(UNAME), Darwin)
../lib/$(SOLIBRARY): $(OBJS)
	@mkdir -p ../lib
	@$(CC) -shared -Wl,-undefined -Wl,dynamic_lookup $(OBJS) -o ../lib/$(SOLIBRARY) $(LDFLAGS)
else
../lib/$(SOLIBRARY): $(OBJS)
	@mkdir -p ../lib
	@$(CC) -shared -o ../lib/$(SOLIBRARY) -Wl,--whole-archive $(OBJS) -Wl,--no-whole-archive $(LDFLAGS)
endif

#
# Generate CTAGS
#

../tags: $(SOURCES) $(INCLUDES)
	@etags -o $@.base $^;\
	cp $@.base $@;\
	for pre in R s c d z; do \
		sed "s/Sprimme/$${pre}primme/g" $@.base >> $@;\
	done

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
deps: $(INCLUDES) $(SOURCES)
	@echo "# This file is generated automatically. Please don't modify" > deps
	@
	@# Add dependencies for automatic generated headers. We ignore headers named after the
	@# source file to avoid cyclic dependencies.
	@tools/dependencies $(INCLUDE) $(SOURCES) --ignore-self-headers -i'include[^\n]+cyclic' -e.h >> deps
	@
	@# Added dependencies on headers for building the objects
	@tools/dependencies $(INCLUDE) $(SOURCES) -e*.o >> deps
	@
	@# Added dependencies on headers from headers not generated automatically
	@tools/dependencies $(INCLUDE) $(INCLUDES) -i'include[^\n]+cyclic' >> deps

include deps

