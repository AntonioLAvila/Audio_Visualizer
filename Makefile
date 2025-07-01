# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -I./src -DSONG_DIR=\"music/snow.wav\"
LDFLAGS := -lsndfile -lportaudio -lsfml-graphics -lsfml-window -lsfml-system -lfftw3

# Directories
SRC_DIR := src
BUILD_DIR := build
BIN := app

# Source and object files
SRCS := $(SRC_DIR)/viz.cpp $(SRC_DIR)/util.cpp
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Default target
all: $(BIN)

# Link final binary
$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compile .cpp to .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(BIN)
