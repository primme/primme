include ../Make_flags

all:
	$(PYTHON) ./setup.py build_ext -i

install:
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

.PHONY : all clean veryclean test install
