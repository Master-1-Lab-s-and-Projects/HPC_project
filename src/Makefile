TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall -O3 -fopenmp#-O0 -g

.PHONY: clean all default cleanall

BUILD_DIR=../build

SRC = main.c datastructure.c exact_cover.c utility.c matrix_parser.c
OBJ_DEPS = $(SRC:.c=.o)
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(OBJ_DEPS))

default: $(TARGET)

all: default

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

$(BUILD_DIR)/%.o: %.c
	bear -a $(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJ)

cleanall: clean
	-rm -f $(TARGET)
