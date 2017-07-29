# Makefile for most code functionality

# References
#   GNU Make manual

# clean: Will remove all files from data and error folder

# C Compiler
CC := gcc

# C standard
CSTD := -std=c11

# Directories for the program
SRCDIR  := src
INCDIR  := include
OBJDIR  := obj
DEPDIR  := depend
EXECDIR := bin

# Name given to the Exectuable generated. Add prefix
# will add the file name
EXECUTABLE := DGSolver.exe
EXECUTABLE := $(addprefix $(EXECDIR)/,$(EXECUTABLE))

LOCAL_INC := -I./include 
LIBS := -framework Accelerate 
INCS := $(LOCAL_INC) 

# Here, wildcard will create a list of all files that have a given
# pattern that is provided. This way we will get all src and header files.

# get the list of all .c files in src
SOURCES := $(wildcard $(SRCDIR)/*.c)

# get the list of all .h files in includes
HEADERS := $(wildcard $(INCDIR)/*.h)

# Here we are doing path replacement. That is, any pattern
# that has the form wiht % sign in sources directory
# has everything replaced other than the percent to 
# the thing on the right.
OBJECTS := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DEPENDS := $(SOURCES:$(SRCDIR)/%.c=$(DEPDIR)/%.d)

# Virtual path to directories that make will search for a prereq file
# that was not found. 
VPATH = ./src:./include

# @ = suppresses the echoing of the command (printing twice)
# $^ = place all prerequisites in the line with spaces in between them

# Compile executable file (Default goal)
$(EXECUTABLE) : $(OBJECTS)
	@echo
	@echo Creating/updating: $@
	@$(CC) -o $@ $^ $(INCS) $(LIBS)
	@echo

# Include dependencies (Must be placed after default goal otherwise some 
# random file will be the default target). Here we are including all the 
# dependency lines by placing the lines in each dependency file here which
# gives us how each object file relates to other header files. Included as the 
# target in each rule is also the dependency file which must be regenerated, like the 
# object file, if any prerequisites change.
include $(DEPENDS)

# include Depends will create a list of rules with the targets as object and
# depend files and prerequisites being the files they both depend on. No recipe has
# been defined yet however. Although a target can have multiple rules, it can only have 
# one recipe. When make executes, all prereqs from all rules are combined and the ONE
# recipe that the target depends on (if a prepreq changes) is run if needed. So, from
# here to the end, specify the recipe to generate the object file and depend file, but
# no prereqs are needed here since they are defined in the depends portion.

# Using patterns in the rule will treat each form with the given stem (%) 
# using the format specified (i.e. its the same as if each file was copied and pasted
# in % and the whole rule was copied and pasted for each file)

# $< = The name of the first prerequisite.
# %@ = Variable for the target

# Create object files
$(OBJDIR)/%.o : %.c
	@echo Creating/updating: $@
	@$(CC) $(CSTD) -c -o $@ $< $(INCS)

# Create the recipes for how to create the dependency file for each source file here. 
# Recall % is a pattern element so each dependency file's matching .c file can 
# be found easily using pattern.

# MM = do not include files found in system files (built in libraries that are linked)
# MG = -MG assumes missing header files are generated files and adds them to the dependency list without raising an error.

$(DEPDIR)/%.d : %.c
	@$(CC) -MM -MG $< > $@; # Use first prerequisite only as other prereqs are included from the existing %.d file
	@sed -i -e 's|.*:|$(OBJDIR)/$*.o $(DEPDIR)/$*.d:|' $@
	@echo Creating/updating: $@

# Phony make functions
.PHONY: clean
clean:
	rm $(EXECUTABLE) $(OBJECTS) $(DEPENDS)




