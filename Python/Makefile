# SWIG
GENERATED  = primme_wrap.cxx primme_wrap.h Primme.py

all: $(GENERATED) primmew.h primmew.cxx Primme.py setup.py
	./setup.py build_ext -i

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

.PHONY : all clean veryclean test swig
