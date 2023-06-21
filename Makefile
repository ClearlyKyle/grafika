# Build settings
CC=gcc

CFLAGS=-g -Wall -mconsole

# SDL options
CC_SDL=$(shell pkg-config --cflags --libs sdl2)

SOURCE = 3d.c
EXEC = 3d

OBJ_TEST_SRC = objtest.c obj.h

all:build

build:
	@echo [BUILD] Building project...
	$(CC) $(CFLAGS) $(SOURCE) -o $(EXEC) $(CC_SDL)
	@echo [BUILD] Building complete!

build_run:build
	@echo [RUN] Building project...
	$(EXEC)
	@echo [RUN] Finished!

run:build
	@echo [RUN] Running...
	./$(EXEC)

clean:
	@echo [CLEAN] Cleaning project...
	-del -fR *.exe *.o
	@echo [CLEAN] Clean completed!

3dbuild:
	@echo [BUILD] Building 3d project...
	$(CC) $(CFLAGS) $(SOURCE) -o $(EXEC) $(CC_SDL)
	@echo [BUILD] Building complete!

objbuild:
	@echo [BUILD] Building obj test...
	$(CC) $(CFLAGS) $(OBJ_TEST_SRC) -o objtest $(CC_SDL)

objtest:objbuild
	@echo [RUN]
	./objtest.exe
	
3d:3dbuild
	@echo [RUN] Running...
	./$(EXEC)