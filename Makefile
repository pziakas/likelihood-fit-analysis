# Makefile for Maximum Likelihood Fit Project using ROOT

CXX = g++
CXXFLAGS = -std=c++17 -O2
ROOTFLAGS = `root-config --cflags --libs` -lMinuit

SRC = src/main.cpp src/likelihood_fit.C
OUT = likelihood_fit

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(ROOTFLAGS) -o $(OUT)

clean:
	rm -f $(OUT)