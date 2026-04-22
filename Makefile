CXXFLAGS+=-std=c++20 -Wall -Wextra -Wconversion -Wshadow
LDLIBS+=-lgsl -lm -lfmt
#LINK.o=${LINK.cc}

all: plots doc.pdf

doc.pdf: doc.tex
	latexmk -xelatex doc.tex

plots: run
	./run
	gnuplot fig6.gnuplot
	gnuplot fig3a.gnuplot


clean:
	rm -f \
		doc.out \
		doc.fls \
		doc.aux \
		doc.xdv \
		doc.log \
		doc.fdb_latexmk \
		*.csv \
		run.o \
		run \

squeaky_clean: clean
	rm -f \
		*.png \
		*.pdf \
