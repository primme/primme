include ../Make_flags

# SWIG
GENERATED  = primme_wrap.cxx primme_wrap.h Primme.py

all: $(GENERATED) primmew.h Primme.py setup.py
	$(PYTHON) setup.py build_ext -i

install: $(GENERATED) primmew.h Primme.py setup.py
	$(PYTHON) setup.py install

# Test target: run the tests
test: examples.py
	$(PYTHON) examples.py

all_tests: tests.py examples.py
	$(PYTHON) tests.py

# Update $(GENERATED) files
swig: primme.i numpy.i wrappers.py primmew.h pyfragments.swg
	swig -c++ -Wextra -python primme.i

# Clean target
clean:
	$(RM) -rf build
	$(RM) -f _Primme*

veryclean: clean
	$(RM) $(GENERATED)

.PHONY : all clean veryclean test swig install
