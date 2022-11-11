rm -rf test_english
make nuclear && make test_english
./test_english test/english.100MB
rm -rf test_english
