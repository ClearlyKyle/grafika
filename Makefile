#
CC := gcc
PKG_CONFIG := pkg-config

CFLAGS := -g -std=gnu11 -ffast-math

# CFLAGS += -gdwarf-2 # DrMemory things

CFLAGS += -fopenmp -march=native
CFLAGS += -mconsole

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
CFLAGS += -Wjump-misses-init
#CFLAGS += -Wcast-align=strict
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

# ignore errors
CFLAGS += -Wno-unused-function

# release flags
CFLAGS_RELEASE := -O3 -DNDEBUG -mconsole -fopenmp -msse4.1

# pkg-config library names
LIB_NAMES := sdl2 SDL2_ttf SDL2_image

LIB_CFLAGS  := $(shell $(PKG_CONFIG) --cflags $(LIB_NAMES))
LIB_LDFLAGS := $(shell $(PKG_CONFIG) --libs-only-L --libs-only-l $(LIB_NAMES)) -lm

INC_DIRS 	:= ./src ./src/graphics/#

CFLAGS         += $(addprefix -I, $(INC_DIRS))
CFLAGS_RELEASE += $(addprefix -I, $(INC_DIRS))

SRC_DIR := src#
BIN_DIR := bin#

MAIN_TARGETS := main terrain

all: $(addprefix $(BIN_DIR)/,$(MAIN_TARGETS:=.exe))

$(BIN_DIR)/%.exe: $(BIN_DIR)/%.o | $(BIN_DIR)
	@echo "Linking $@..."
	@$(CC) -fopenmp $(LDFLAGS) $^ -o $@ $(LIB_LDFLAGS)
	@echo "Successfully built $@"

$(BIN_DIR)/%.o: $(SRC_DIR)/%.c | $(BIN_DIR)
	@echo "Compiling main $<..."
	@mkdir -p $(dir $@)
	@$(CC) $(CFLAGS) $(LIB_CFLAGS) -MMD -MP -c $< -o $@

# include auto-generated dependencies for MAIN_TARGETS only
-include $(addprefix $(BIN_DIR)/,$(MAIN_TARGETS:=.d))

$(MAIN_TARGETS): %: $(BIN_DIR)/%.exe

$(BIN_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BIN_DIR)
	@echo "Cleaned build directory"

flags:
	@echo $(CFLAGS)

.PHONY: all clean $(MAIN_TARGETS)