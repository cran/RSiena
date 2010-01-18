## This should work on Win64 but it not currently built.
include $(R_HOME)/etc$(R_ARCH)/Makeconf

%.o: %.rc
	$(RESCOMP) -i $< -o $@

all: siena.exe sienaC.exe

install: all
	cp siena.exe sienaC.exe ../../inst

siena.exe: siena.o siena_rc.o
	$(CC) -s -mwindows -o $@ siena.o siena_rc.o

sienaC.exe: siena.o siena_rc.o
	$(CC) -s -o $@ siena.o siena_rc.o

clean:
	@rm -f *.exe *.o ../../inst/*.exe
