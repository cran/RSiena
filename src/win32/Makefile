%.o: %.rc
	windres -i $< -o $@

all: siena.exe sienaC.exe

install: all
	cp siena.exe sienaC.exe ../../inst

siena.exe: siena.o siena_rc.o
	gcc -s -mwindows -o $@ siena.o siena_rc.o

sienaC.exe: siena.o siena_rc.o
	gcc -s -o $@ siena.o siena_rc.o

clean:
	@rm -f *.exe *.o ../../inst/*.exe
