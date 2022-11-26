rm -rf test_english
make nuclear && make test_english
./test_english test/kernel
rm -rf test_english
