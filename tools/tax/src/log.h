#pragma once
#include <iostream>
#include <ctime>
#include <iomanip>
#include <assert.h>

#define LOG(expr) { \
    auto t = time(nullptr); \
    auto lt = std::localtime(&t); \
    char timestamp[24]; \
    auto written = strftime(timestamp, 24, "%Y-%m-%d %H:%M:%S", lt); \
    assert(written < 24); \
    std::cerr << timestamp << "\t" << expr << std::endl; \
}
