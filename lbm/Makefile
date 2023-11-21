defaults: all

include ../config.mk

#CPPFLAGS = -I.
#CXXFLAGS = -O3 -fopenmp -funroll-loops -std=c++11
#LDFLAGS = -lgomp -lrt -lpng

include ../flags.mk

CPP_SOURCES = $(wildcard *.cpp)
CU_SOURCES = $(wildcard *.cu)
#SOURCES = $(CPP_SOURCES) $(CU_SOURCES)
SIMS = $(wildcard sim?*.cu)
TARGETS = $(SIMS:%.cu=%)

# OBJS are *.o except sim?*.o
OBJS = $(filter-out $(SIMS:%.cu=%.o), $(CU_SOURCES:%.cu=%.o))  $(CPP_SOURCES:%.cpp=%.o)


all: $(TARGETS)

clean:
	rm -f *.o *.d $(TARGETS)

%.o: %.cpp
	$(CXX) -MD -MP $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

%.o: %.cu
ifeq ($(use_CUDA),yes)
	$(NVCC) -M $(NVCCFLAGS) -o $(@:%.o=%.d) $<
	$(NVCC) $(CPPFLAGS) $(NVCCFLAGS) -c -o $@ $<
else
#use CXX instead of NVCC --- important for using lgomp and fopenmp for CPU based LBM
	$(CXX) -MD -MP $(CPPFLAGS) $(CXXFLAGS) -xc++ -c -o $@ $<
endif

$(TARGETS): % : %.o $(OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $< $(OBJS) $(LDFLAGS)

# All sources depend on config.mk
#%.cpp: ../config.mk Makefile

-include $(CPP_SOURCES:%.cpp=%.d)
-include $(CU_SOURCES:%.cu=%.d)
