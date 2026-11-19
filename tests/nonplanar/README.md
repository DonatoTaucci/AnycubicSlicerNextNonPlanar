Standalone tests for `src/libslic3r/GCode/NonPlanarCore.*` (no other dependency):

    g++ -std=c++17 -O1 -fsanitize=address,undefined -I../../src/libslic3r/GCode \
        test_nonplanar_core.cpp ../../src/libslic3r/GCode/NonPlanarCore.cpp -o t && ./t
