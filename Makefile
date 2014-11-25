gimli: greedy.c
	gcc -Wall  -o $@ greedy.c -lm

gimli_optimized : greedy.c
	gcc -Wall  -o $@ greedy.c -lm -O3
	

gimli_profile : greedy.c
	gcc -Wall  -o $@ greedy.c -lm -g -pg


out.gimli.2: gimli 
	./gimli G199.sample > out.gimli
	awk '$$NF==2' out.gimli > out.gimli.2

test: out.gimli.2.ref out.gimli.2 
	diff $^

clean:
	rm -f gimli gimli_profile gimli_optimized out.gimli.2

