
app := nbody

sources := $(wildcard *.c)
headers := $(wildcard *.h)
src := $(sources) $(headers)
OBJECTS = $(sources:.c=.o)
OFLAGS = -O3
STD = -std=c99

all: $(app)

CFLAGS := -c -g
OMPFLAG := -fopenmp

$(app).o: $(src)
	$(CC) $(sources) $(STD) $(CFLAGS) $(OMPFLAG)

$(app): $(app).o
	$(CC) $(OBJECTS) $(STD) -o $@ -lm

.PHONY: clean, tmpfiles

o3:
	$(CC) $(sources) $(STD) $(OFLAGS) $(OMPFLAG) -lm -o $(app) 

clean:
	@rm *.o $(app) 2> /dev/null || true

tmpfiles:
	@rm Population* stage* 2> /dev/null || true
