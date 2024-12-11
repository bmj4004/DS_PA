# Makefile


main : main.o
	gcc -o $@ $^ -lGL -lGLU -lglut -lpng -lm
	rm *.o

main.o : main.c
	gcc -c $^ -o $@

clean :
	rm main
