include ../Make_flags

all: update
	$(PYTHON) ./setup.py build_ext -i

install:
	$(PYTHON) ./setup.py install

update: primme.c
	rm -rf src
	mkdir -p src/primme src/primme/include src/primme/linalg src/primme/eigs src/primme/svds
	cd ..; cp include/*.h Python/src/primme/include
	cd ../src; find . -name '*.h' -exec cp {} ../Python/src/primme/{} \;
	cd ../src; find . -name '*.c' -exec cp {} ../Python/src/primme/{} \;

# Test target: run the tests
test: examples.py
	$(PYTHON) examples.py

all_tests: tests.py examples.py
	$(PYTHON) tests.py

primme.c: primme.pyx
	cython primme.pyx

# Clean target
clean:
	$(RM) -rf build primme.*so primme.egg-info __pycache__ src dist

veryclean: clean

.PHONY : all clean veryclean test install
