# Make sure there are no spaces after setting variables
# Check for correct indentation!

# Build settings
CC = gcc

CFLAGS := -std=gnu11 -g -Wall -mconsole -fopenmp

EXEC 		:= ModelViewer#
OUTPUT_DIR 	:= bin#
SRC_DIR 	:= src#

INC_DIRS 	:= src deps#

# pkg-config Library names
LIB_NAMES := sdl2

# Library flags and linking
CFLAGS += $(shell pkg-config --cflags $(LIB_NAMES))
LDFLAGS := $(shell pkg-config --libs $(LIB_NAMES))

# Add the include directories
CFLAGS += $(addprefix -I, $(INC_DIRS))

# List of source files, all of the .c files
SOURCES := $(wildcard $(SRC_DIR)/*.c)

# Generate object file names based on source file names
# Convert source files to object file names
# put the % part of the first section, into the % part of the second section
# using the data from the last section
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(OUTPUT_DIR)/%.o,$(SOURCES))

# Default target
all: $(OUTPUT_DIR)/$(EXEC)

# Rule to build the executable
$(OUTPUT_DIR)/$(EXEC): $(OBJECTS)
	@echo Building Exe : $^
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Build the object files
$(OBJECTS): $(SOURCES)
	@echo Building Oject : $<
	$(CC) $(CFLAGS) -c $< -o $@


debug_vars:
	@echo $(OUTPUT_DIR)/$(EXEC)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo $(CFLAGS)
	@echo $(LDFLAGS)


# without phony, make assumes "build" is a file that needs to be made
.PHONY: build
build: $(OUTPUT_DIR)/$(EXEC)

.PHONY: run
run: build
	$(OUTPUT_DIR)/$(EXEC)

.PHONY: clean
clean:
	rm -f $(OUTPUT_DIR)/*.exe $(OUTPUT_DIR)/*.o
	@echo [CLEAN] Clean completed!