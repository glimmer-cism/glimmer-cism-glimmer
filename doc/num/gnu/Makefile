PLOTS=w_profile.pslatex wt_sigma.pslatex

EPS =  

all:			$(PLOTS) $(EPS)

.gp.eps:				
			gnuplot $<
			ps2epsi $@ .tmp
			mv .tmp  $@

.gp.pslatex:		
			gnuplot $<
			sed 's#psfile=#psfile=\\dir/gnu/#' $@ > .tmp
			mv .tmp  $@

.SUFFIXES:		.gp .pslatex .eps

clean:
			rm -f $(PLOTS) $(EPS)
