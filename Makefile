
giml: giml.c
	gcc -Wall  -o $@ $^ -lm

giml_optimized : giml.c
	gcc -Wall  -o $@ $^ -lm -O3

giml_profile : giml.c
	gcc -Wall  -o $@ $^ -lm -g -pg

giml_static: giml.c
	gcc -Wall  -o $@ $^ -lm -static	


.PHONY : clean 

clean:
	rm -f giml giml_profile giml_static giml_optimized 
