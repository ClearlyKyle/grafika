# Make sure there are no spaces after setting variables
# Check for correct indentation!

# Build settings
CC = C:\msys64\mingw64\bin\gcc.exe
PKG_CONFIG = C:\msys64\mingw64\bin\pkg-config.exe

CFLAGS := -g -std=gnu11 -m64

# CFLAGS += -gdwarf-2 # DrMemory things

CFLAGS += -mconsole -fopenmp -msse4.1

CFLAGS += -Wall -Wextra -Wpedantic
#CFLAGS += -Werror # treat warnings as errors

## Security warnings
CFLAGS += -Wformat=2
CFLAGS += -Wformat-overflow=2
CFLAGS += -Wformat-truncation=2
CFLAGS += -Wformat-security
CFLAGS += -Wnull-dereference
CFLAGS += -Wstack-protector
CFLAGS += -Wstrict-overflow=3
CFLAGS += -Wtrampolines
CFLAGS += -Walloca # -Walloca-larger-than=1048576
CFLAGS += -Wvla # -Wvla-larger-than=1048576
CFLAGS += -Warray-bounds=2
CFLAGS += -Wshift-overflow=2
CFLAGS += -Wcast-qual
CFLAGS += -Wstringop-overflow=4
CFLAGS += -Wconversion
#CFLAGS += -Wtraditional-conversion
#CFLAGS += -Warith-conversion

CFLAGS += -Wlogical-op
CFLAGS += -Wduplicated-cond
CFLAGS += -Wduplicated-branches

## GCC 12
CFLAGS += -ftrivial-auto-var-init=zero
CFLAGS += -Wtrivial-auto-var-init
CFLAGS += -Wuse-after-free=3

## Extra flags
CFLAGS += -Wformat-signedness
CFLAGS += -Wshadow
CFLAGS += -Wstrict-overflow=4
CFLAGS += -Wswitch-default
CFLAGS += -Wswitch-enum
CFLAGS += -Wcast-align=strict
CFLAGS += -Wjump-misses-init
#CFLAGS += -Wstack-usage=<byte-size> # Warn if stack usage might exceed <byte-size>
#CFLAGS += -Wstrict-prototypes
#CFLAGS += -Wundef

## Compilation flags
CFLAGS += -fstack-protector-strong
CFLAGS += -fstack-clash-protection
CFLAGS += -fPIE
CFLAGS += -fcf-protection=full # <full|return|branch>

## Glibc flags
#CFLAGS += -D_FORTIFY_SOURCE=2 # GCC < 12, will enable additional security features of the GNU libc when calling memory and string handling functions Ref
CFLAGS += -D_FORTIFY_SOURCE=3 # GCC 12,  will try to detect overflows in variable length variables

#CFLAGS += -fsanitize=address \
#-fsanitize=pointer-compare \pkg
#-fsanitize=pointer-subtract\
#-fno-omit-frame-pointer \
#-fsanitize=undefined \
#-fsanitize=bounds-strict \
#-fsanitize=float-divide-by-zero \
#-fsanitize=leak\
#-fsanitize=float-cast-overflow
#export ASAN_OPTIONS=strict_string_checks=1:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1:detect_invalid_pointer_pairs=2

# release flags
CFLAGS_RELEASE := -O3 -DNDEBUG -mconsole -fopenmp -msse4.1

EXEC 		:= ModelViewer#
OUTPUT_DIR 	:= bin#
SRC_DIR 	:= src#

INC_DIRS 	:= ./src ./deps ./src/graphics/#

# pkg-config Library names
LIB_NAMES := sdl2 sdl2_ttf sdl2_image

# Library flags and linking
LIB_CFLAGS  := $(shell $(PKG_CONFIG) --cflags $(LIB_NAMES)) -mconsole
LIB_LDFLAGS := $(shell $(PKG_CONFIG) --libs-only-L --libs-only-l $(LIB_NAMES))

# Add the include directories
CFLAGS         += $(addprefix -I, $(INC_DIRS))
CFLAGS_RELEASE += $(addprefix -I, $(INC_DIRS))

# List of source files, all of the .c files
SOURCES := $(wildcard $(SRC_DIR)/*.c)

# Generate object file names based on source file names
# Convert source files to object file names
# put the % part of the first section, into the % part of the second section
# using the data from the last section
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(OUTPUT_DIR)/%.o,$(SOURCES))
DEPENDS := $(patsubst $(SRC_DIR)/%.c,$(OUTPUT_DIR)/%.d,$(SOURCES))

# Default target
all: $(OUTPUT_DIR)/$(EXEC)

# Rule to build the executable
$(OUTPUT_DIR)/$(EXEC): $(OBJECTS)
	@echo Building Exe : $^
	$(CC) $(CFLAGS) $^ -o $@ $(LIB_LDFLAGS)

# Build the object files
$(OBJECTS): $(SOURCES)
	@echo Building Oject : $<
	$(CC) $(CFLAGS) $(LIB_CFLAGS) -MMD -MP -c $< -o $@

# without phony, make assumes "build" is a file that needs to be made
.PHONY: build
build: $(OUTPUT_DIR)/$(EXEC)

# Release target
.PHONY: release
release:
	@make build CFLAGS="$(CFLAGS_RELEASE)"

.PHONY: run
run:
	$(OUTPUT_DIR)/$(EXEC)

.PHONY: clean
clean:
	rm -f $(OUTPUT_DIR)/*.exe $(OUTPUT_DIR)/*.o $(OUTPUT_DIR)/*.d
	@echo [CLEAN] Clean completed!

-include $(DEPENDS)

debug_vars:
	@echo $(OUTPUT_DIR)/$(EXEC)
	@echo $(SOURCES)
	@echo $(OBJECTS)
	@echo $(CFLAGS)
	@echo $(LIB_LDFLAGS)
	@echo $(LIB_CFLAGS)
