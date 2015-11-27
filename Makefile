SOURCES=src/EC3NP.cpp src/ECLNP.cpp src/EC2NP.cpp src/smallx.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: executable

$(OBJECTS): %.o : %.cpp

.cpp.o:
	g++ -Wall -O3 -c $< -o $@

executable: $(OBJECTS)
	g++ -Wall -pedantic -Wextra -O3 $^ src/main.cpp -o bin/hpl-eval -lcln -lginac -g

clean:
	rm -f $(OBJECTS) bin/hpl-eval
