make jacobi
make gs
./jacobi > jacobi.txt
./gs > gs.txt

diff jacobi.txt gs.txt