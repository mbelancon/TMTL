# make            -> tllam_fast, tmhamtl_fast e libtlcore (.so/.dll)
# Windows (MinGW): make EXE=.exe LIB=tlcore.dll
CC ?= gcc
CFLAGS ?= -O3 -std=c99 -D_DEFAULT_SOURCE -fopenmp
EXE ?=
LIB ?= libtlcore.so

all: tllam_fast$(EXE) tmhamtl_fast$(EXE) $(LIB)

tllam_fast$(EXE): tllam_fast.c tlcore.c tlcore.h tlcli_common.h
	$(CC) $(CFLAGS) tllam_fast.c tlcore.c -o $@ -lm

tmhamtl_fast$(EXE): tmhamtl_fast.c tlcore.c tlcore.h tlcli_common.h
	$(CC) $(CFLAGS) tmhamtl_fast.c tlcore.c -o $@ -lm

$(LIB): tlcore.c tlcore.h
	$(CC) $(CFLAGS) -shared -fPIC tlcore.c -o $@ -lm

clean:
	rm -f tllam_fast$(EXE) tmhamtl_fast$(EXE) $(LIB)
