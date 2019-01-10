
RM := rm -rf

CXX = mpicxx
CXXFLAGS = -std=c++11

BUILDDIR := $(shell mkdir -p build; echo build)

SRCS += \
src/main.cpp \
src/TreeNode.cpp \
src/Utilities.cpp 

OBJS += \
$(BUILDDIR)/main.o \
$(BUILDDIR)/TreeNode.o \
$(BUILDDIR)/Utilities.o

$(BUILDDIR)/%.o: src/%.cpp
	@echo 'Building file: $<'
	$(CXX) $(CXXFLAGS) $< -c -o $@
	@echo 'Finished building: $<'
	@echo ' '

all: Barnes-Hut

Barnes-Hut: $(OBJS)
	@echo 'Building target: $@'
	$(CXX) $(CXXFLAGS) $(OBJS) -o "Barnes-Hut"
	@echo 'Finished building target: $@'
	@echo ' '
	
clean:
	-$(RM) $(OBJS)
	-$(RM) $(BUILDDIR)
	-$(RM) Barnes-Hut
	-@echo ' '
	