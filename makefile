all:  docs install

docs:
	R -e "devtools::document()"
vig:
	R -e "devtools::build_vignettes()"

build:
	R CMD build --no-build-vignettes .

install:
	R CMD INSTALL --no-test-load .


test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

