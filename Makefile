# Define the build directory
BUILD_DIR = build
EXECUTABLE = LIE-ALG

# Phony targets are not real files
.PHONY: all prepare configure build run clean

# Default target
all: configure build run

# Prepare the build directory
prepare:
	rm -r $(BUILD_DIR)
	mkdir $(BUILD_DIR)

# Configure the project
configure:
	cd $(BUILD_DIR) && cmake ..

# Build the project
build:
	cd $(BUILD_DIR) && cmake --build .

# Run the executable
run:
	cd $(BUILD_DIR) && ./$(EXECUTABLE)

# Clean the build directory
clean:
	rm -r $(BUILD_DIR)

exp:
	pandoc -s output.md -o output.pdf
