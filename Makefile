CC ?= cc

CPPFLAGS := -Iinclude
CFLAGS ?= -O2
CFLAGS += -std=c11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow

BUILD_DIR := build
LIB_OBJECT := $(BUILD_DIR)/affine_align.o
CLI_OBJECT := $(BUILD_DIR)/main.o
TEST_OBJECT := $(BUILD_DIR)/test_affine_align.o
CLI := $(BUILD_DIR)/affine-align
TEST_BINARY := $(BUILD_DIR)/affine-align-tests

.PHONY: all test sanitize clean

all: $(CLI)

test: $(TEST_BINARY) $(CLI)
	$(TEST_BINARY)
	$(CLI) --fasta examples/pair.fasta --format json >/dev/null

sanitize:
	$(MAKE) BUILD_DIR=build-sanitize \
		CFLAGS="-O1 -g -std=c11 -Wall -Wextra -Wpedantic -Wconversion \
		-Wshadow -fsanitize=address,undefined -fno-omit-frame-pointer" \
		LDFLAGS="-fsanitize=address,undefined" test

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(LIB_OBJECT): src/affine_align.c include/affine_align.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(CLI_OBJECT): src/main.c include/affine_align.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TEST_OBJECT): tests/test_affine_align.c include/affine_align.h | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(CLI): $(LIB_OBJECT) $(CLI_OBJECT)
	$(CC) $(LDFLAGS) $^ -o $@

$(TEST_BINARY): $(LIB_OBJECT) $(TEST_OBJECT)
	$(CC) $(LDFLAGS) $^ -o $@

clean:
	$(RM) -r $(BUILD_DIR)
