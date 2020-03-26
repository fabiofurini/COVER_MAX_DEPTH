	
	
	# COMPILER

CXX= g++

# FILES

OBJS =	src/Main.o src/global_functions_KP.o src/global_functions_COVER.o src/DP.o 


# CPLEX VERSION (LIBS and INCLUDE files)

#cplex 12.7 HERE

CPLEXLIBDIR =  /home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/lib/x86-64_linux/static_pic
LP_INCLUDE= /home/fabio/ILOG/CPLEX_Studio_AcademicResearch129/cplex/include/ilcplex

#CPLEXLIBDIR =  /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/lib/x86-64_osx/static_pic
#LP_INCLUDE= /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/include/ilcplex


#CPLEXLIBDIR =  /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/lib/x86-64_osx/static_pic
#LP_INCLUDE = /Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/include/ilcplex

#CPLEXLIBDIR = /opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic
#LP_INCLUDE = /opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ilcplex


# Nothing should be changed

LP_LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -lm -lpthread -ldl

DBG= -O3
#DBG= -g


#DEFS = $(OS_VERSION) $(COMPILER) $(LP_SOLVERS)

INCDIR = -I. -I$(LP_INCLUDE)

#COMPILER FLAGS

CXXFLAGS =  $(DBG) $(INCDIR) 

.c.o:
	gcc -c $(CXXFLAGS) $< -o $@

LDLIBS = $(LP_LIBS)

all:COVER

COVER: $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDLIBS) 

$(OBJS): Makefile

clean:
	rm -f $(OBJS) rm COVER

	

