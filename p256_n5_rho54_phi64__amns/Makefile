OPTIONS= -O3 -Wall -g -lgmp -lcrypto

CC= gcc-10

OPTIONS+=$(MACROS)

main : main.o useful_functs.o
	$(CC) -o main main.o useful_functs.o $(OPTIONS)


useful_functs.o : useful_functs.c useful_functs.h 
	$(CC) -c useful_functs.c $(OPTIONS)


main.o : useful_functs.o
	$(CC) -c main.c $(OPTIONS)


	
all : main


clean :
	rm *.o -fv *~

cleanall : clean
	rm main -f

