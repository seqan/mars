#include <gtest/gtest.h>

#include <seqan3/std/iterator>
//#include <seqan3/std/ranges>
//#include <string_view>
#include <vector>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "motif.hpp"
#include "multiple_alignment.hpp"

TEST(Motif, Detection)
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

    std::vector<mars::StemloopMotif> motifs = mars::detect_stemloops(bpseq, plevel);
    EXPECT_EQ(motifs.size(), 2);
    EXPECT_EQ(motifs[0].bounds, (mars::Coordinate{27, 47}));
    EXPECT_EQ(motifs[1].bounds, (mars::Coordinate{54, 68}));
    EXPECT_EQ(motifs[0].uid, 0);
    EXPECT_EQ(motifs[1].uid, 1);
}

TEST(Motif, AnalyzeStemLoop)
{
    std::vector<int> bpseq{76,75,74,73,72,71,70,-1,-1,-1,
                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                           -1,-1,-1,-1,-1,-1,-1,47,46,45,
                           44,43,-1,-1,-1,-1,-1,-1,-1,-1,
                           -1,-1,-1,31,30,29,28,27,-1,-1,
                           -1,-1,-1,-1,68,67,66,65,-1,-1,
                           -1,-1,-1,-1,-1,57,56,55,54,-1,
                           6, 5, 4, 3, 2, 1, 0,-1};

    mars::Msa msa{};
    msa.sequences.resize(5);
    using std::ranges::copy;
    copy(std::string_view{"gcuuuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguagaa-aacu--a-ua---cguuuuaaagcu"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(msa.sequences[0]));
    copy(std::string_view{"acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaaau-uuuggugcaacuccaaauaaaagua"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(msa.sequences[1]));
    copy(std::string_view{"gcgguuguggcgaag-ugguuaacgcaccagauuguggcucuggcacuc----guggguucgauucccaucaaucgcc"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(msa.sequences[2]));
    copy(std::string_view{"gcucauguagc-ucaguugguagagcacacccu----ugguaagggugaggucagcgguucaaauccgcucaugagcu"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(msa.sequences[3]));
    copy(std::string_view{"guuucuguagu-ugaau---uacaacgaugauu----uuucaugucauuggu-cgcaguugaaugcuguguagaaaua"}
         | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(msa.sequences[4]));

    mars::StemloopMotif motif{0, {27, 47}};
    motif.analyze(msa, bpseq);
    EXPECT_EQ(motif.bounds, (mars::Coordinate{27ul, 47ul}));
    EXPECT_EQ(motif.length.min, 17);
    EXPECT_EQ(motif.length.max, 21);
    EXPECT_FLOAT_EQ(motif.length.mean, 17.8);
    EXPECT_EQ(motif.elements.size(), 2);

    // check the stem
    EXPECT_TRUE(std::holds_alternative<mars::StemElement>(motif.elements[0]));
    mars::StemElement const & stem = std::get<mars::StemElement>(motif.elements[0]);
    EXPECT_EQ(stem.length.min, 10);
    EXPECT_EQ(stem.length.max, 10);
    EXPECT_FLOAT_EQ(stem.length.mean, 10);
    EXPECT_EQ(stem.profile.size(), 5);
    EXPECT_RANGE_EQ(stem.profile[0].quantities(), (std::array<float, 16>{0,0,0,2,0,0,1,1,0,0,0,0,1,0,0,0}));
    EXPECT_RANGE_EQ(stem.profile[1].quantities(), (std::array<float, 16>{0,0,0,1,0,1,1,0,0,0,0,0,2,0,0,0}));
    EXPECT_EQ(stem.gaps.size(), 5);
    for (auto const & map : stem.gaps)
    {
        EXPECT_TRUE(map.empty());
    }

    // check the loop
    EXPECT_TRUE(std::holds_alternative<mars::LoopElement>(motif.elements[1]));
    mars::LoopElement loop = std::get<mars::LoopElement>(motif.elements[1]);
    EXPECT_FALSE(loop.is_5prime);
    EXPECT_EQ(loop.length.min, 7);
    EXPECT_EQ(loop.length.max, 11);
    EXPECT_FLOAT_EQ(loop.length.mean, 7.8);
    EXPECT_EQ(loop.profile.size(), 11);
    EXPECT_RANGE_EQ(loop.profile[10].quantities(), (std::array<float, 4>{0,2,0,3}));
    EXPECT_RANGE_EQ(loop.profile[9].quantities(), (std::array<float, 4>{0,0,0,1}));
    EXPECT_EQ(loop.gaps.size(), 11);
    for (int idx : {0, 1, 2, 3, 4, 5, 6, 7, 8, 10})
    {
        EXPECT_TRUE(loop.gaps[idx].empty());
    }
    EXPECT_EQ(loop.gaps[9].size(), 1);
    auto const & gap_entry = loop.gaps[9].begin();
    EXPECT_EQ(gap_entry->first, 4);
    EXPECT_EQ(gap_entry->second, 4);
}
