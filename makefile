make: compile run

compile:
	gcc -Wall main.c -o main -O3 -lm
run:
	./main
limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)
