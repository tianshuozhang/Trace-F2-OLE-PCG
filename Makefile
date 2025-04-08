PCG_TARGET = ./bin/pcg
DPF_TARGET = ./bin/dpf

CC = gcc
CXX = g++ 
CFLAGS += -std=c99 -O3 -I./include -I./libs/fft/include -I./libs/gf64/include -I./libs/gf128/include -I./libs/tri-dpf/include -I/usr/include/openssl/ -I/usr/local/include/emp-ot/ -I/usr/local/include/emp-tool/
CXXFLAGS += $(CFLAGS) -msse4.1 
LDFLAGS = -march=native -lcrypto -lssl -lm -maes -ffast-math -lpthread

# 公共源文件（不包含main文件）
FFT_SRC = $(filter-out ./libs/fft/src/test.c, $(wildcard ./libs/fft/src/*.c))
DPF_SRC = $(filter-out ./libs/tri-dpf/src/test.c, $(wildcard ./libs/tri-dpf/src/*.c))
GF128_SRC = $(filter-out ./libs/gf128/src/test.c, $(wildcard ./libs/gf128/src/*.c))
GF64_SRC = $(filter-out ./libs/gf64/src/test.c, $(wildcard ./libs/gf64/src/*.c))
OTDPF_CPP_SRC = ./libs/tri-dpf/src/otdpf.cpp


# PCG专用源文件
PCG_C_SRC = $(filter-out ./src/main.c, $(wildcard ./src/*.c))
PCG_CPP_SRC = ./src/main.cpp

# DPF专用源文件
DPF_CPP_SRC = $(wildcard ./test/*.cpp)


# 公共对象文件
COMMON_OBJS = $(FFT_SRC:.c=.o) \
              $(DPF_SRC:.c=.o) \
              $(GF128_SRC:.c=.o) \
              $(GF64_SRC:.c=.o) \
              $(OTDPF_CPP_SRC:.cpp=.o)

# PCG专用对象文件
PCG_OBJS = $(PCG_C_SRC:.c=.o) \
           $(PCG_CPP_SRC:.cpp=.o)




# DPF专用对象文件
DPF_OBJS = $(DPF_CPP_SRC:.cpp=.o)


all: $(PCG_TARGET) $(DPF_TARGET) 

$(PCG_TARGET): $(COMMON_OBJS) $(PCG_OBJS)
	@mkdir -p ./bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DPF_TARGET): $(COMMON_OBJS) $(DPF_OBJS)
	@mkdir -p ./bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)




# C编译规则
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# C++编译规则
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f \
		libs/fft/src/*.o \
		libs/tri-dpf/src/*.o \
		libs/gf128/src/*.o \
		libs/gf64/src/*.o \
		src/*.o \
		test/*.o \
		$(PCG_TARGET) \
		$(DPF_TARGET)

.PHONY: all clean