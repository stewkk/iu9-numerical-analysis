.PHONY: build

all: build

build:
	@latexmk -f -pdf -output-directory=build -shell-escape ./report.tex
	@latexmk -f -pdf -output-directory=build -shell-escape ./report2.tex
	@latexmk -f -pdf -output-directory=build -shell-escape ./report3.tex
	@latexmk -f -pdf -output-directory=build -shell-escape ./report4.tex
	@latexmk -f -pdf -output-directory=build -shell-escape ./report5.tex
	@latexmk -f -pdf -output-directory=build -shell-escape ./report6.tex

build-presentation:
	@latexmk -f -pdf -output-directory=build -shell-escape ./presentation.tex

clean:
	@rm -rf build

build-with-biblio: build
	@cp ./biblio.bib ./build
	@cd build && biber report

view:
	@evince ./build/report.pdf

view-presentation:
	@evince ./build/presentation.pdf
