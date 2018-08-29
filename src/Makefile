CC	= gcc
CFLAGS	= -O3
LDFLAGS	= 
INCLUDES = 
LIBS	= -lm
TARGET	= fastMatch
OBJS	= fastGMatch.o calFeature.o utilities.o

# Commands
all:	$(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS) core

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

calFeature.o: fastGMatch.h calFeature.h utilities.h
fastGMatch.o: fastGMatch.h calFeature.h utilities.h
utilities.o: fastGMatch.h calFeature.h utilities.h
