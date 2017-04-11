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


update:
	cp -r ../src/eigs ../src/include ../src/linalg ../src/svds src/primme
	cp -r ../include inst
	R -e "require(Rcpp); compileAttributes();"
	R -e "require(devtools); document()"
	R -e 'require(rmarkdown); render("README.Rmd",output_format="all")'

doc:
	find . -name '*.pdf' -exec rm {} \;
	R -e "require(devtools); document()"
	find . -name '*.pdf' -exec rm {} \;
	R CMD Rd2pdf .
	find . -name '*.pdf' -exec rm {} \;
	
install: update
	R CMD INSTALL .

clean:
	find . -name '*.pdf' -exec rm {} \;
	find . -name '*.o' -exec rm {} \;
	find . -name '*.so' -exec rm {} \;
	find . -name '*.orig' -exec rm {} \;
	rm -f src/primme/libprimme.a

.PHONY: update doc install clean
