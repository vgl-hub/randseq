CXX = g++
INCLUDE_DIR = -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS)

TARGET = mytool
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o
LIBS = -lz
LDFLAGS = -pthread

MODULE1_DIR := $(CURDIR)/submodule1

OBJS := main
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(MODULE1_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

all: head

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@
	
.PHONY: module1
module1:
	$(MAKE) -j -C $(MODULE1_DIR) CXXFLAGS="$(CXXFLAGS)"
	
.PHONY: module2
module2:
	$(MAKE) -j -C $(MODULE2_DIR) CXXFLAGS="$(CXXFLAGS)"
	
$(BUILD):
	-mkdir -p $@
	
clean:
	$(MAKE) -j -C $(SUBMODULE1_SUBDIR) clean
	$(MAKE) -j -C $(SUBMODULE2_SUBDIR) clean
	$(RM) -r build
