CPP	= g++
CFLAGS	= -O3
LDFLAGS	= 
INCLUDES = 
LIBS	= -lm
TARGET	= fastGMatch
OBJS	= fastGMatch.o extFeature.o utilities.o

# Commands
all:	$(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS) core

.cpp.o:
	$(CPP) $(CFLAGS) $(INCLUDES) -c $<

calFeature.o: fastGMatch.h extFeature.h utilities.h parameter.h
fastGMatch.o: fastGMatch.h extFeature.h utilities.h parameter.h
utilities.o: fastGMatch.h extFeature.h utilities.h parameter.h
