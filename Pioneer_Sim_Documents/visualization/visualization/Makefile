CC = g++
CFLAGS = -Wall
PROG = flight

SRCS = ac3d.c texture.c render.c Imageloader.cpp flight_mae_154A.cpp

ifeq ($(OS),Windows_NT)
    LIBS = -lfreeglut -lglew32 -lopengl32 -lglu32
else
    LIBS = -framework OpenGL -framework GLUT
endif

all: $(PROG)

$(PROG): $(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)