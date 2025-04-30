# Makefile for lab1.cpp

# 编译器
CXX = g++

# 编译选项
CXXFLAGS = -Wall -g

# 源文件
SRC = lab1.cpp

# 目标可执行文件
TARGET = lab1

# 默认目标
all: $(TARGET)

# 规则：如何编译源文件
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# 清理目标
clean:
	rm -f $(TARGET) *.png *.dot *.txt
