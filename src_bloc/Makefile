TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall -O3

.PHONY: clean all default cleanall

BUILD_DIR=../build/build_bloc

SRC = exact_cover.c 
OBJ_DEPS = $(SRC:.c=.o)
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(OBJ_DEPS))

ifneq (, $(shell which bear 2>/dev/null))
BEAR_CMD = bear -a
endif

default: $(TARGET)

all: default

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(@D)
	$(BEAR_CMD) $(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJ)

cleanall: clean
	-rm -f $(TARGET)
