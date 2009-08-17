#/usr/bin/bash

.PHONY: all
all : bd gfrd hardbody irreversible mapk pushpull rebind reversible tf

DATA=data/run/files.pvd data/run/files/*

.PHONY: bd
bd:
	echo
	cd bd_test && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python irr.py data/output.out 0.000001 10

.PHONY: cylinders
cylinders:
	echo
	cd cylinders && rm -rf $(DATA) && LOGLEVEL=DEBUG PYTHONPATH=$(PYTHONPATH):../../ python run.py

.PHONY: gfrd
gfrd:
	echo
	# cd gfrd && rm -rf $(DATA) && LOGLEVEL=ERROR  
	# PYTHONPATH=$PYTHONPATH:../../ python single-gfrd.py

.PHONY: hardbody
hardbody:
	echo
	cd hardbody && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run_single.py 1 1 100

.PHONY: irreversible
irreversible:
	echo
	cd irreversible && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py data/output.out 1 100

.PHONY: mapk
mapk:
	echo
	cd mapk && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python mapk.py

.PHONY: pushpull
pushpull:
	echo
	cd pushpull && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python pushpull.py 0.1 0.1 10 20 30 1 'single' 0.00001

.PHONY: rebind
rebind:
	echo
	cd rebind && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py 1 20 0 1

.PHONY: reversible
reversible:
	echo
	cd reversible && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py data/output.out 0.05 5


.PHONY: single
single:
	echo
	cd single && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py single.0.out 5e-4 3e-8 100


.PHONY: tf
tf:
	echo
	cd tf && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python tf.py

