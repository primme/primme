# ------------------------------------------------------------------------------
# The next action may help developers to maintain the R package:
#
#   make update: copy PRIMME source files in include and src to this package;
#                also update automatic generated code by Rcpp and roxygen2.
#
#        doc: generate and show the PDF documentation.
#
#        clean: remove object files.
#
#        install: install the local version of the package.


update: clean src/primme/zhegv.c
	rm -rf inst/include src/primme/include src/primme/linalg src/primme/eigs src/primme/svds
	mkdir -p src/primme inst/include inst/include/c_interface src/primme/include src/primme/linalg src/primme/eigs src/primme/svds
	cd ..; cp include/*.h R/inst/include/c_interface
	cd ../src; find . -name '*.h' -exec cp {} ../R/src/primme/{} \;
	cd ../src; find . -name '*.c' -exec cp {} ../R/src/primme/{}pp \;
	find src/primme -name '*.cpp' -exec sed -i -e '/#define THIS_FILE/ s/\.c"/.cpp"/' {} \;
	R -e "require(Rcpp); compileAttributes();"
	touch inst/include/PRIMME_types.h
	R -e "require(devtools); load_all(); document()"
	#R -e 'require(rmarkdown); render("README.Rmd",output_format="all")'
	make clean

src/primme/zhegv.c:
	mkdir -p src/primme
	# R LAPACK may not have zhegv; R installations often use a complete LAPACK, but this is required to pass the r-cran package test
	wget -q "http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=lapack%2Fcomplex16%2Fzhegv.f" -O - | ./f2clapack.sh > $@

doc:
	find . -name '*.pdf' -exec rm {} \;
	R -e "require(devtools); load_all(); document()"
	find . -name '*.pdf' -exec rm {} \;
	R CMD Rd2pdf .
	find . -name '*.pdf' -exec rm {} \;
	
install: update
	R CMD INSTALL .

check: update
	R -e "require(devtools); load_all(); check()"

clean:
	find . -name '*.pdf' -exec rm {} \;
	find . -name '*.o' -exec rm {} \;
	find . -name '*.a' -exec rm {} \;
	find . -name '*.so' -exec rm {} \;
	find . -name '*.orig' -exec rm {} \;
	./cleanup

.PHONY: update doc install clean check
