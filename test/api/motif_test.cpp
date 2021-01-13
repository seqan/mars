#include <gtest/gtest.h>

#include <vector>

#include <seqan3/test/expect_range_eq.hpp>

#include "motif.hpp"

TEST(Stemloop, Detection)
{
    std::vector<int> bpseq{76,75,74,73,72,71,70,-1,-1,-1,
                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                           -1,-1,-1,-1,-1,-1,-1,47,46,45,
                           44,43,-1,-1,-1,-1,-1,-1,-1,-1,
                           -1,-1,-1,31,30,29,28,27,-1,-1,
                           -1,-1,-1,-1,68,67,66,65,-1,-1,
                           -1,-1,-1,-1,-1,57,56,55,54,-1,
                            6, 5, 4, 3, 2, 1, 0,-1};
    std::vector<int> plevel{ 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1, 0, 0, 0,
                             0, 0,-1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1, 0, 0, 0, 0, 0,-1,-1,
                            -1,-1,-1,-1, 0, 0, 0, 0,-1,-1,
                            -1,-1,-1,-1,-1, 0, 0, 0, 0,-1,
                             0, 0, 0, 0, 0, 0, 0,-1};

    mars::stemloop_type stemloops = mars::detect_stem_loops(bpseq, plevel);

    mars::stemloop_type expected {{27,47}, {54,68}};
    EXPECT_RANGE_EQ(stemloops, expected);
}
