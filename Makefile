
gimli: gimli.c
	gcc -Wall  -o $@ $^ -lm

gimli_optimized : gimli.c
	gcc -Wall  -o $@ $^ -lm -O3

gimli_profile : gimli.c
	gcc -Wall  -o $@ $^ -lm -g -pg

gimli_static: greedy.c
	gcc -Wall  -o $@ $^ -lm -static	


.PHONY : clean 

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized 
