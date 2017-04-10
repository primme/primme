# SWIG
GENERATED  = primme_wrap.cxx primme_wrap.h Primme.py

all: $(GENERATED) primmew.h primmew.cxx Primme.py setup.py
	./setup.py build_ext -i

install: $(GENERATED) primmew.h primmew.cxx Primme.py setup.py
	./setup.py install

# Test target: run the tests
test: tests.py examples.py
	python tests.py

# Update $(GENERATED) files
swig: primme.i numpy.i wrappers.py primmew.h pyfragments.swg
	swig -c++ -Wextra -python primme.i

# Clean target
clean:
	$(RM) -rf build
	$(RM) -f _Primme*

veryclean: clean
	$(RM) $(GENERATED)

# Generate the directory primme with the source code of PRIMME
update:
	rm -rf primme
	mkdir -p primme primme/include primme/src primme/src/include primme/src/linalg primme/src/eigs primme/src/svds
	cd ..; find include src -name '*.h' -exec cp {} Python/primme/{} \;
	for f in `cd ../src; find . -name '*.c'` ; do \
		cp ../src/$$f primme/src/$$f; \
		echo "#define USE_DOUBLE\n#define F77UNDERSCORE\n#include \"`basename $$f`\"" > primme/src/$${f%%.c}double.c; \
		echo "#define USE_DOUBLECOMPLEX\n#define F77UNDERSCORE\n#include \"`basename $$f`\"" > primme/src/$${f%%.c}doublecomplex.c; \
		echo "#define USE_FLOAT\n#define F77UNDERSCORE\n#include \"`basename $$f`\"" > primme/src/$${f%%.c}float.c; \
		echo "#define USE_FLOATCOMPLEX\n#define F77UNDERSCORE\n#include \"`basename $$f`\"" > primme/src/$${f%%.c}floatcomplex.c; \
	done

.PHONY : all clean veryclean test swig install update
