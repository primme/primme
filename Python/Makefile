include ../Make_flags

all:
	./setup.py build_ext -i

install:
	./setup.py install

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
	cython --cplus primme.pyx

# Generate the directory primme with the source code of PRIMME
update: primme.cpp
	rm -rf primme
	mkdir -p primme primme/include primme/src primme/src/include primme/src/linalg primme/src/eigs primme/src/svds
	cd ..; find include src -name '*.h' -exec cp {} Python/primme/{} \;
	cd ..; find include src -name '*.c' -exec cp {} Python/primme/{}pp \;
	find primme -name '*.cpp' -exec sed -i -e '/#define THIS_FILE/ s/\.c"/.cpp"/' {} \;

.PHONY : all clean veryclean test install update
