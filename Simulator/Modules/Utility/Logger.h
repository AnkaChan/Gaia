#pragma once

#define DEBUG_LVL_WARNING 0 
#define DEBUG_LVL_INFO    1
#define DEBUG_LVL_DEBUG   2
#define DEBUG_LVL_DEBUG_VEBOSE   3
#define DEBUG_LVL_DEBUG_VEBOSE_2   4

namespace GAIA {
    inline void debugInfoGen(const int debugLvlSet, int debugLvl, std::function<void()> ops, bool trigger = true) {
        if (debugLvl <= debugLvlSet && trigger)
        {
            ops();
        }
    };

    //debugInfoGen(curPhysics, 3, [&]() {
    //    std::cout << "Debug Info." << std::endl;
    //    });
}