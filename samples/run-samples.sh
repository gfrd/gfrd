#!/bin/bash


cd bd_test
echo
echo 'bd_test'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python irr.py data/output.out 0.000001 1
# Todo: sometimes I get this: 
#
# Traceback (most recent call last):
#   File "irr.py", line 62, in <module>
#     run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
#   File "irr.py", line 10, in run
#     d, t = singlerun( T )
#   File "irr.py", line 55, in singlerun
#     distance = s.distance( particleB.getPos(), particleA.getPos() )
#   File "/home/thomas/amolf/epdp-0.3/gfrdbase.py", line 239, in getPos
#     return pool.positions[ pool.indexMap[ self.serial ] ]
# KeyError: 0


cd ../gfrd
echo
echo 'gfrd doesnt work. Todo.'
#LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python single-gfrd.py


cd ../hardbody
echo
echo 'hardbody'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python run_single.py 1 1 100


cd ../irreversible
echo
echo 'irreversible'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python run.py data/output.out 1 100


cd ../mapk
echo
echo 'mapk'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python mapk.py


cd ../pushpull
echo 'pushpull'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python pushpull.py 0.1 0.1 10 20 30 1 'single' 0.00001


cd ../rebind
echo
echo 'rebind'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python run.py 1 20 0 1


cd ../reversible
echo
echo 'reversible'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python run.py data/output.out 0.05 5


cd ../tf
echo
echo 'tf'
LOGLEVEL=ERROR PYTHONPATH=$PYTHONPATH:../../ python tf.py


echo
cd ../






### Graveyard

#Dirlist=$(find -mindepth 1 -maxdepth 1 -type d)
#for direc in $Dirlist ;
#do
#    cd $direc
#    echo $direc
#    make
#    echo cd ../
#done

