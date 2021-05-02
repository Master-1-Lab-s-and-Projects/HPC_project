TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall -O3

.PHONY: clean all default cleanall

BUILD_DIR=../build_dynamic

SRC = main.c datastructure.c exact_cover.c utility.c matrix_parser.c
OBJ_DEPS = $(SRC:.c=.o)
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(OBJ_DEPS))

default: $(TARGET)

all: default

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJ)

cleanall: clean
	-rm -f $(TARGET)
