rm -rf test_english
rm -rf test/kernel
cd test
wget http://pizzachili.dcc.uchile.cl/repcorpus/real/kernel.gz
gzip -d kernel.gz
rm -rf kernel.gz
cd ..
make nuclear && make test_english
./test_english test/kernel
rm -rf test_english
rm -rf test/kernel