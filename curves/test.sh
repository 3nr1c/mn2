if [ $(grep -i TODO curves.c | wc -l) -gt 0 ]
then
	echo -e "\n\t\x1B[93;5mWARNING: pending TODOs in curves.c\x1B[0m\n"
fi

make curves
./curves
