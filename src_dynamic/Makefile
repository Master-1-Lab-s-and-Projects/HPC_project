TARGET = exact_cover.exe
CC = mpicc
CFLAGS = -Wall -DDEBUG

BUILD_DIR=../build/build_$(subst src_,,$(notdir $(shell pwd)))

SRC = main.c datastructure.c exact_cover.c utility.c matrix_parser.c
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

.PHONY: clean all default cleanall
