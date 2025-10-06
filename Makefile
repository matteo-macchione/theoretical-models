# ===========================
# Makefile per theoreticalsnew.cpp con CBL
# ===========================

# Compilatore
CXX = g++
CXXFLAGS = -O2 -std=c++11 -Wall -Wextra

# Percorsi include
INCLUDES = \
    -I/home/matteo/CosmoBolognaLib-master/Headers \
    -I/home/matteo/CosmoBolognaLib-master/External/Eigen/eigen-3.4.0

# Percorsi librerie
LIBPATH = -L/home/matteo/CosmoBolognaLib-master

# Librerie da linkare
LIBS = -lCOSM -lFUNC -lFUNCGRID -lWRAP_LIB -lgsl -lgslcblas -lm

# Nome eseguibile
TARGET = theoreticalsnew

# File sorgente
SRC = theoreticalsnew.cpp

# ===========================
# Regole Makefile
# ===========================

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $(TARGET) $(LIBPATH) $(LIBS)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean



