
# This Makefile does not really 'make' anything.

# Running 'make' will run all sample scripts. Used for debugging.


.PHONY: all
all : bd dimer example hardbody irreversible mapk pushpull rebind reversible single tf

DATA=data/run/files.pvd data/run/static.pvd data/run/files/
PYTHONPATH_NEW=${PYTHONPATH}:../../

.PHONY: bd
bd:
	echo
	cd bd_test && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python irr.py data/output.out 0.000001 10

.PHONY: dimer
dimer:
	echo
	cd dimer && echo 'dimer is old'

.PHONY: example
example:
	echo
	cd example && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=$(PYTHONPATH_NEW) python example.py

.PHONY: hardbody
hardbody:
	echo
	cd hardbody && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python run_single.py 1 1 100

.PHONY: irreversible
irreversible:
	echo
	cd irreversible && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python run.py data/output.out 1 100

.PHONY: mapk
mapk:
	echo
	cd mapk && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python mapk.py

.PHONY: pushpull
pushpull:
	echo
	cd pushpull && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python pushpull.py 0.1 0.1 10 20 30 1 'single' 0.00001

.PHONY: rebind
rebind:
	echo
	cd rebind && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python run.py 1 20 0 1

.PHONY: reversible
reversible:
	echo
	cd reversible && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python run.py data/output.out 0.05 5


.PHONY: single
single:
	echo
	cd single && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python run.py single.0.out 5e-4 3e-8 100


.PHONY: tf
tf:
	echo
	cd tf && rm -rf $(DATA) && LOGLEVEL=ERROR PYTHONPATH=${PYTHONPATH_NEW} python tf.py

