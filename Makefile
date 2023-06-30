# Build settings
CC=gcc

CFLAGS=-g -Wall -mconsole

# pkg-config Library names
LIB_NAMES := sdl2

# Library flags and linking
CFLAGS += $(shell pkg-config --cflags $(LIB_NAMES))
LDFLAGS := $(shell pkg-config --libs $(LIB_NAMES))

# Executable name
EXEC := ModelViewer

INCLUDE_FOLDER = src deps
CFLAGS += $(addprefix -I, $(INCLUDE_FOLDER))

# Specify the output directory
OUTPUT_DIR := bin

# Source directory
SRC_DIR := src

# List of source files
SRCS := $(wildcard $(SRC_DIR)/*.c $(SRC_DIR)/*.h)

# Generate object file names based on source file names
OBJS := $(patsubst $(SRC_DIR)/%.c,$(OUTPUT_DIR)/%.o,$(SRCS))

# Default target
all: $(OUTPUT_DIR)/$(EXEC)

# Rule to build the executable
$(OUTPUT_DIR)/$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Rule to build object files
$(OUTPUT_DIR)/%.o: $(SRC_DIR)/%.c | $(OUTPUT_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Phony target to run the executable
run_main: build_main clean
	$(OUTPUT_DIR)/$(EXEC)

.PHONY: build_main
build_main: clean $(OUTPUT_DIR)/$(EXEC)
	@echo [BUILD] Building complete!

.PHONY: clean
clean:
	rm -f $(OUTPUT_DIR)/*.exe $(OUTPUT_DIR)/*.o
	@echo [CLEAN] Clean completed!