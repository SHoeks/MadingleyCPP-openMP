# Source directories separated by space
# Example ./ src1/ src2/
SRCDIR=./ Data/ Input/ Model/ Output/ Tools/
# Directory where object files will be placed
OBJDIR=../build/
# Include directories separated by space
# Example: include1/ include2/
#INCDIR=./ Data/ Input/ Model/ Output/ Tools/  /usr/local/netcdf-cxx4-4.3.0/include /usr/local/netcdf-4.4.1/include 
INCDIR=./ Data/ Input/ Model/ Output/ Tools/  
# Libraries separated by space
#LIBS= $(SUBLIBS)  -Wl,-rpath,$(/usr/local/netcdf-cxx4-4.3.0/lib)  -L/usr/local/netcdf-cxx4-4.3.0/lib -L/usr/local/netcdf-4.4.1/lib -L/usr/lib -lnetcdf_c++4 -lnetcdf -lhdf5_hl_cpp -lhdf5_cpp -lz -ldl -lpthread 
LIBS= $(SUBLIBS)  -Wl,-rpath, -L/usr/lib -lnetcdf_c++4 -lnetcdf -lhdf5_hl_cpp -lhdf5_cpp -lz -ldl -lpthread 
# Directory where binary file will be placed
BINDIR=../dist/
# Name of the result file
TARGET=madingley
# Compiler
CXX=g++

# Retrive list of the source files
SRC=$(wildcard $(addsuffix *.cpp,$(SRCDIR)))
# Generate list of the object files
OBJ=$(addprefix $(OBJDIR), $(patsubst %.cpp, %.o, $(notdir $(SRC))))

VPATH=$(SRCDIR)

# Compilation flags
CXXFLAGS=-O2 -std=c++11 -fopenmp

$(TARGET) : $(OBJ)
	@echo Linking...
	@mkdir -p $(BINDIR)
	@$(CXX) $(CXXFLAGS) -o $(BINDIR)$@ $(OBJ) $(LIBS) -fopenmp
	@cp -R -n Model_setup/input $(BINDIR)
	@echo Madingley compilation complete.

$(OBJDIR)%.o : %.cpp
	@echo Compiling $< in $@...
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) $(addprefix -I,$(INCDIR)) -c -o $@ $< -fopenmp

clean :
	@$(RM) -r $(OBJDIR)
	@$(RM) -r $(BINDIR)
