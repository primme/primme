
all:
	./setup.py build_ext -i

install:
	./setup.py install

# Test target: run the tests
test: tests.py examples.py
	python tests.py

# Clean target
clean:
	$(RM) -rf build primme.so

veryclean: clean

.PHONY : all clean veryclean test install
