
RM      = rm -f
CC     ?= gcc
CFLAGS ?= -Wall -Wextra

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c,%.o,$(SRC))


default: su

clean:
	$(RM) su $(OBJ)

.PHONY: clean


su: $(OBJ)
	$(CC) -o $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

