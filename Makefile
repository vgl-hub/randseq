CXX = g++
INCLUDE_DIR = -I./include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = mytool
BUILD = build/bin
SOURCE = src
INCLUDE = include
LDFLAGS :=

SUBMODULE1_SUBDIR := $(CURDIR)/submodule1
SUBMODULE1_LIBSFILES := $(GFASTATS_SUBDIR)/$(SOURCE)/* $(SUBMODULE1_SUBDIR)/$(INCLUDE)/*

SUBMODULE2_SUBDIR := $(CURDIR)/submodule2
SUBMODULE2_LIBSFILES := $(GFALIGN_SUBDIR)/$(SOURCE)/* $(SUBMODULE2_SUBDIR)/$(INCLUDE)/*

main: $(SOURCE)/main.cpp $(SUBMODULE1_LIBSFILES) $(SUBMODULE2_LIBSFILES) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOURCE)/main.cpp -o $(BUILD)/$(TARGET)

$(SUBMODULE1_LIBSFILES): submodule1
	@# Do nothing
	
$(SUBMODULE2_LIBSFILES): submodule2
	@# Do nothing

.PHONY: submodule1
submodule1:
	$(MAKE) -j -C $(SUBMODULE1_SUBDIR)
	
.PHONY: submodule2
submodule2:
	$(MAKE) -j -C $(SUBMODULE2_SUBDIR)
	
$(BUILD):
	-mkdir -p $@
	
clean:
	$(MAKE) -j -C $(SUBMODULE1_SUBDIR) clean
	$(MAKE) -j -C $(SUBMODULE2_SUBDIR) clean
	$(RM) -r build
