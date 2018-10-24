function removeTest {
	sed 's|/\*\* BEGIN TEST \*\*/|\n&|g;s|/\*\* END TEST \*\*/|&\n|g' $1 | sed '/\/\*/,/*\//d' > $2
}

if [ ! -d "build" ]; then
	mkdir build
fi;

# syntax: /** BEGIN TEST **/ (...) /** END TEST **/
removeTest linear.c build/linear.c