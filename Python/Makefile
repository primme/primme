include ../Make_flags

all: primme.cpp primme-gpu.cpp
	$(PYTHON) ./setup.py build_ext -i

install: primme.cpp primme-gpu.cpp
	$(PYTHON) ./setup.py install

# Test target: run the tests
test: examples.py
	$(PYTHON) examples.py

all_tests: tests.py examples.py
	$(PYTHON) tests.py

# Clean target
clean:
	$(RM) -rf build primme.*so

veryclean: clean

primme.cpp: primme.pyx
	cython --cplus -E BUILD_WITH_MAGMA=False -o primme.cpp primme.pyx

primme-gpu.cpp: primme.pyx
	cython --cplus -E BUILD_WITH_MAGMA=True -o primme-gpu.cpp primme.pyx

.PHONY : all clean veryclean test install
