make linear

./linear

shasum -a 256 b.txt | diff test-assets/b_sha256.txt -
diff jacobi.txt gs.txt