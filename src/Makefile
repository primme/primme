
include ../Make_flags

VPATH = ../include:include:linalg:eigs:svds

find_sources = $(wildcard $(dir)/*.c)
DIRS := linalg eigs svds
SOURCES := $(foreach dir,$(DIRS),$(find_sources))

ARCHS := float floatcomplex double doublecomplex
OBJS := $(foreach arch,$(ARCHS),$(patsubst %.c,%$(arch).o,$(SOURCES)))

INCLUDE += -I../include -Iinclude

MACRO_HEADERS := \
   ../include/primme.h \
   include/template.h \

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
   eigs/locking.h \
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

%double.o: %.c %.d
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_DOUBLE -c $< -o $@

%doublecomplex.o: %.c %.d
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_DOUBLECOMPLEX -c $< -o $@

%float.o: %.c %.d
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_FLOAT -c $< -o $@

%floatcomplex.o: %.c %.d
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDE) -DUSE_FLOATCOMPLEX -c $< -o $@

#
# Archive object files in the full library
#

../lib/$(LIBRARY): $(OBJS)
	@mkdir -p ../lib;\
	ar rv $@ $(OBJS);\
	$(RANLIB) $@

../lib/$(SOLIBRARY): ../lib/$(LIBRARY)
	$(CC) -shared -o ../lib/$(SOLIBRARY) -Wl,--whole-archive ../lib/$(LIBRARY) -Wl,--no-whole-archive

#
# Generate CTAGS
#

../tags: $(SOURCES) $(MACRO_HEADERS)
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
.PHONY: auto_headers clean %.d

#
# Dependencies
#

include/template.h: primme.h
include/numerical.h: template.h blaslapack.h auxiliary.h
linalg/blaslapack.d: blaslapack.h template.h blaslapack_private.h
linalg/auxiliary.d: auxiliary.h template.h blaslapack.h
linalg/wtime.d: wtime.h

eigs/auxiliary_eigs.d: auxiliary.h const.h numerical.h globalsum.h wtime.h
eigs/convergence.d: convergence.h const.h numerical.h ortho.h auxiliary_eigs.h
eigs/correction.d: correction.h const.h numerical.h inner_solve.h globalsum.h auxiliary_eigs.h
eigs/factorize.d: factorize.h numerical.h
eigs/globalsum.d: globalsum.h numerical.h
eigs/init.d: init.h numerical.h update_projection.h update_W.h ortho.h factorize.h wtime.h auxiliary_eigs.h
eigs/inner_solve.d: inner_solve.h numerical.h inner_solve.h factorize.h update_W.h globalsum.h wtime.h auxiliary_eigs.h
eigs/locking.d: locking.h const.h numerical.h convergence.h auxiliary_eigs.h restart.h 
eigs/main_iter.d: main_iter.h const.h wtime.h numerical.h main_iter_private.h convergence.h correction.h init.h ortho.h restart.h solve_projection.h update_projection.h update_W.h globalsum.h auxiliary_eigs.h
eigs/ortho.d: ortho.h numerical.h globalsum.h const.h
eigs/primme.d: const.h wtime.h numerical.h convergence.h correction.h init.h ortho.h restart.h solve_projection.h update_projection.h update_W.h primme_interface.h
eigs/primme_f77.d: primme_f77_private.h primme_interface.h notemplate.h
eigs/primme_f77_private.h: template.h
eigs/primme_interface.d: template.h const.h primme_interface.h notemplate.h
eigs/restart.d: restart.h const.h numerical.h locking.h ortho.h solve_projection.h factorize.h update_projection.h update_W.h convergence.h globalsum.h auxiliary_eigs.h
eigs/solve_projection.d: solve_projection.h const.h numerical.h ortho.h
eigs/update_projection.d: update_projection.h const.h numerical.h globalsum.h
eigs/update_W.d: update_W.h numerical.h ortho.h auxiliary_eigs.h wtime.h

svds/primme_svds.d: numerical.h wtime.h primme_svds_interface.h primme_interface.h
svds/primme_svds_f77.d: primme_svds_f77_private.h primme_svds_interface.h notemplate.h
svds/primme_svds_f77_private.h: template.h
svds/primme_svds_interface.d: numerical.h primme_interface.h primme_svds_interface.h notemplate.h

