TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall 

.PHONY: clean all default cleanall

SRC = main.c datastructure.c exact_cover.c utility.c matrix_parser.c
OBJ = $(SRC:.c=.o)

default: $(TARGET)

all: default

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

%.o: %.c
	bear -a $(CC) -c -o $@ $<

clean:
	-rm -f $(OBJ)

cleanall: clean
	-rm -f $(TARGET)
