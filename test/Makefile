TARGET = ./bin/pcg
CC = gcc
CXX = g++
CFLAGS += -std=c99 -O3 -I ../include -I../libs/tri-dpf/include -I/usr/include/openssl/ -I ../libs/emp-ot
CXXFLAGS += $(CFLAGS)
LDFLAGS = -march=native -lcrypto -lssl -lm -maes -ffast-math

# 源文件定义
DPF_SRC = $(filter-out ../libs/tri-dpf/src/test.c, $(wildcard ../libs/tri-dpf/src/*.c))
PCG_CPP_SRC = ./main.cpp
MAIN_CPP_SRC = ../libs/tri-dpf/src/otdpf.cpp

OBJECTS = $(MAIN_CPP_SRC:.cpp=.o) \
          $(DPF_SRC:.c=.o) \
          $(PCG_CPP_SRC:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p ./bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)  # 使用g++链接

# C编译规则
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# C++编译规则
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f ../libs/tri-dpf/src/*.o *.o $(TARGET) $(OBJECTS)

.PHONY: all clean