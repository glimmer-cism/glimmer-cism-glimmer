PSTEX = scale.pstex_t grid.pstex_t
EPS   = 

all:		$(PSTEX) $(EPS)

.fig.pstex_t:	
		fig2dev -Lpstex_t -p $*.pstex -m 1.000000 $< $@
		#sed 's#file=#file=\\dir/figs/#' $@ > .tmp
		sed 's#includegraphics{#includegraphics{\\dir/figs/#' $@ > .tmp
		mv .tmp  $@
		fig2dev -Lpstex -m 1.000000 -n $*.pstex $< $*.pstex

.fig.eps:
		fig2dev -Lps -m 1.000000 -n $*.eps $< $*.eps

.SUFFIXES:	.fig .pstex .pstex_t .eps

clean:
		rm *.pstex *.pstex_t *.eps
