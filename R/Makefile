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
	mkdir -p src/primme inst/include src/primme/include src/primme/linalg src/primme/eigs src/primme/svds
	cd ..; find include -name '*.h' -exec cp {} R/inst/{} \;
	cd ../src; find . -name '*.[ch]' -exec cp {} ../R/src/primme/{} \;
	R -e "require(Rcpp); compileAttributes();"
	R -e "require(devtools); load_all(); document()"
	#R -e 'require(rmarkdown); render("README.Rmd",output_format="all")'
	make clean

src/primme/zhegv.c:
	mkdir -p src/primme
	# R LAPACK may not have zhegv; R installations often use a complete LAPACK, but this is required to pass the r-cran package test
	wget -q "http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=lapack%2Fcomplex16%2Fzhegv.f" -O - | ./f2clapack.sh > src/primme/zhegv.c

doc:
	find . -name '*.pdf' -exec rm {} \;
	R -e "require(devtools); document()"
	find . -name '*.pdf' -exec rm {} \;
	R CMD Rd2pdf .
	find . -name '*.pdf' -exec rm {} \;
	
install: update
	R CMD INSTALL .

check: update
	R -e "require(devtools); check()"

clean:
	find . -name '*.pdf' -exec rm {} \;
	find . -name '*.o' -exec rm {} \;
	find . -name '*.a' -exec rm {} \;
	find . -name '*.so' -exec rm {} \;
	find . -name '*.orig' -exec rm {} \;
	./cleanup

.PHONY: update doc install clean check
