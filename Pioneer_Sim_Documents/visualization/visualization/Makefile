CC = g++
CFLAGS = -Wall
PROG = flight

#SRCS = ac3d.c texture.c render.c vec3f.cpp Imageloader.cpp flight_carrier_test.cpp
SRCS = ac3d.c texture.c render.c Imageloader.cpp flight_mae_154A.cpp
#SRCS = ac3d.c texture.c render.c Imageloader.cpp flight_klz.cpp

ifeq ($(shell uname),Darwin)
	LIBS = -framework OpenGL -framework GLUT
else
	LIBS = -lglut
endif

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)
