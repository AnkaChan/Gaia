#pragma once
#include <chrono>
#include <ctime>

#define TICK(name)\
auto name##_t1 = std::chrono::high_resolution_clock::now();

#define TOCK(name)\
auto name##_t2 = std::chrono::high_resolution_clock::now();\
name = std::chrono:: duration_cast<std::chrono::microseconds>(name##_t2 - name##_t1).count() / 1000.0;

#define TOCK_STRUCT(statisticsStruct, name)\
auto name##_t2 = std::chrono::high_resolution_clock::now();\
statisticsStruct.name += std::chrono:: duration_cast<std::chrono::microseconds>(name##_t2 - name##_t1).count() / 1000.0;

