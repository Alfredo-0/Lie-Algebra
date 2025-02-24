# Define the build directory
BUILD_DIR = build

# OS-specific settings
ifeq ($(OS),Windows_NT)
    RM = rmdir /S /Q
    MKDIR = mkdir
    BUILD_CONFIG = Debug
    EXECUTABLE = LIE-ALG.exe
    EXEC_PREFIX =
else
    RM = rm -r
    MKDIR = mkdir -p
    BUILD_CONFIG =
    EXECUTABLE = LIE-ALG
    EXEC_PREFIX = ./
endif

# Phony targets are not real files
.PHONY: all prepare configure build run clean exp

# Default target
all: configure build run

# Prepare the build directory
prepare:
	-$(RM) $(BUILD_DIR)
	$(MKDIR) $(BUILD_DIR)

# Configure the project
configure:
	cd $(BUILD_DIR) && cmake ..

# Build the project
build:
	cd $(BUILD_DIR) && cmake --build .

# Run the executable
run:
	cd $(BUILD_DIR)/$(BUILD_CONFIG) && $(EXEC_PREFIX)$(EXECUTABLE)

# Clean the build directory
clean:
	-$(RM) $(BUILD_DIR)

# Additional target for converting markdown to PDF
exp:
	pandoc -s output.md -o output.pdf
