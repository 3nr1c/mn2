echo "- Compiling linear.c"
make linear

echo "- Testing 10"
./linear 10
diff jacobi.txt test-assets/10_test.txt
diff gs.txt test-assets/10_test.txt

echo "- Testing 10e6"
./linear
shasum -a 256 b.txt | diff test-assets/b_sha256.txt -
diff jacobi.txt gs.txt