#/usr/bin/bash

.PHONY: all
all : bd gfrd hardbody irreversible mapk pushpull rebind reversible tf



.PHONY: bd
bd:
	echo
	cd bd_test && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python irr.py data/output.out 0.000001 1

.PHONY: cylinders
cylinders:
	echo
	cd cylinders && rm -rf data/run && LOGLEVEL=DEBUG PYTHONPATH=$(PYTHONPATH):../../ python run.py

.PHONY: gfrd
gfrd:
	echo
	#cd gfrd && rm -rf data/* && LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python single-gfrd.py

.PHONY: hardbody
hardbody:
	echo
	cd hardbody && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run_single.py 1 1 100

.PHONY: irreversible
irreversible:
	echo
	cd irreversible && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py data/output.out 1 100

.PHONY: mapk
mapk:
	echo
	cd mapk && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python mapk.py

.PHONY: pushpull
pushpull:
	echo
	cd pushpull && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python pushpull.py 0.1 0.1 10 20 30 1 'single' 0.00001

.PHONY: rebind
rebind:
	echo
	cd rebind && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py 1 20 0 1

.PHONY: reversible
reversible:
	echo
	cd reversible && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python run.py data/output.out 0.05 5

.PHONY: tf
tf:
	echo
	cd tf && rm -rf data/run && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH}:../../ python tf.py







### Graveyard

#Dirlist=$(find -mindepth 1 -maxdepth 1 -type d)
#for direc in $Dirlist ;
#do
#    cd $direc
#    echo $direc
#    make
#    echo cd ../
#done

