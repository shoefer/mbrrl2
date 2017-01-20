BASE   = ../..
NAME   = $(shell basename $(BASE))
OUTPUT = $(BASE)/lib/$(NAME).so
CXXFLAGS += -fPIC

SYS = $(shell uname)
ifneq ($(SYS), Darwin)
LIBS   += -lrt
endif

OPTIM  = fast_debug
#OPTIM  = debug

LAPACK = 0

CXXFLAGS += -DMT_IMPLEMENT_TEMPLATES #causes all *.cpp to implement the used templates

OBJS  = $(shell find . -maxdepth 1 \( -name '*.cpp' -and -not -name '*_t.cpp' -and -not -name '*_instantiate.cpp' \) )
OBJS := $(OBJS:%.cpp=%.o)
OBJS += ../MT/array.o \
	../MT/util.o

include $(BASE)/make-generic
