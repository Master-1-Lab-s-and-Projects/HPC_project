TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall -O0 -g

.PHONY: clean all default cleanall

SRC = main.c datastructure.c exact_cover.c utility.c matrix_parser.c
OBJ = $(SRC:.c=.o)

default: $(TARGET)

all: default

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

%.o: %.c
	bear -a $(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJ)

cleanall: clean
	-rm -f $(TARGET)
