get md5 from
https://launchpad.net/madgraph5
gunzip+untar

./bin/mg5
mg5>install Delphes

Alternatively (note some more recent versions fail due to a missing Pythia.h):

wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.0.9.tar.gz
tar xzf Delphes-3.0.9.tar.gz
cd Delphes-3.0.9
./configure
make
