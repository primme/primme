# SWIG
INTERFACES = primme.i
WRAPPERS   = primme_wrap.cxx primme_wrap.h
PROXIES    = Primme.py

# Default target: build the tests
all: $(WRAPPERS) primmew.h primmew.cxx
	./setup.py build_ext -i

# Test target: run the tests
test: all
	python tests.py

# Rule: %.i -> %_wrap.cxx
%_wrap.cxx %_wrap.h: %.i numpy.i wrappers.py primmew.h
	swig -c++ -Wextra -python $<

# Clean target
clean:
	$(RM) -rf build
	$(RM) -f *.so  $(WRAPPERS) $(PROXIES)*

.PHONY : all clean test
