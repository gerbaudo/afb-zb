SAMPLES = \
	sm_zb \
	sm_zbbar \
	bm_zb \
	bm_zbbar

all : plotHistos

setup :
	@echo "todo: build delphes library"
	@echo "todo: retrieve or link data files"
	@echo "for now link existing ones"
	data/create_links.sh

fillHistos : # todo: add dependece on input ntuple files
	echo $(SAMPLES) | xargs -n 1 -P 4 python/plot_zb.py 
#	the command above runs the four samples in parallel
# 	python/plot_zb.py sm_zb
# 	python/plot_zb.py sm_zbbar
# 	python/plot_zb.py bm_zb
# 	python/plot_zb.py bm_zbbar

plotHistos : fillHistos # todo: depend on output root files instead
	python/plot_zb.py

