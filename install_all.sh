#!/bin/bash
#
## First create a virtual environment:
virtualenv -p python2.7 .venv
#
## Then activate it
. .venv/bin/activate
#
## Then run this script
# ./install_all.sh
#
#
## FUSED-Wake should now be accessible from your new virtual environment
#

echo 'INSTALL OPENMDAO'

pip install numpy scipy

curl http://openmdao.org/releases/0.10.3.2/go-openmdao-0.10.3.2.py | python2

echo 'INSTALL FUSED-Wake'
python setup.py install

