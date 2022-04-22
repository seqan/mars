#include <gtest/gtest.h>

#include <seqan3/std/iterator>
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
    EXPECT_EQ(motifs.size(), 3u);
    EXPECT_EQ(motifs[0].bounds, (mars::Coordinate{27, 47}));
    EXPECT_EQ(motifs[1].bounds, (mars::Coordinate{54, 68}));
    EXPECT_EQ(motifs[2].bounds, (mars::Coordinate{0, 26}));
    EXPECT_EQ(motifs[0].uid, 0);
    EXPECT_EQ(motifs[1].uid, 1);
    EXPECT_EQ(motifs[2].uid, 2);
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
    msa.structure.first = std::move(bpseq);
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
    motif.analyze(msa);
    EXPECT_EQ(motif.bounds, (mars::Coordinate{27ul, 47ul}));
    EXPECT_EQ(motif.length.min, 17);
    EXPECT_EQ(motif.length.max, 21);
    EXPECT_FLOAT_EQ(motif.length.mean, 17.8);
    EXPECT_EQ(motif.elements.size(), 2u);

    // check the stem
    EXPECT_TRUE(std::holds_alternative<mars::StemElement>(motif.elements[1]));
    mars::StemElement const & stem = std::get<mars::StemElement>(motif.elements[1]);
    EXPECT_EQ(stem.prio.size(), 5u);

    // check stem.prio[4] values: {0,0,0,2,0,0,1,1,0,0,0,0,1,0,0,0}
    EXPECT_EQ(stem.prio[4].size(), 4ul);
    auto stemIt = stem.prio[4].cbegin();
    EXPECT_EQ(stemIt->second.to_chars(), (std::pair<char, char>{'C', 'G'}));
    EXPECT_FLOAT_EQ(stemIt->first, log2f(17328 * (1*600+1) * 2. / 600 / 5 / 9913));
    ++stemIt;
    EXPECT_EQ(stemIt->second.to_chars(), (std::pair<char, char>{'U', 'A'}));
    EXPECT_FLOAT_EQ(stemIt->first, log2f(17328 * (1*600+1) * 2. / 600 / 5 / 3975));
    ++stemIt;
    EXPECT_EQ(stemIt->second.to_chars(), (std::pair<char, char>{'A', 'U'}));
    EXPECT_FLOAT_EQ(stemIt->first, log2f(17328 * (2*600+1) * 2. / 600 / 5 / 3975));
    ++stemIt;
    EXPECT_EQ(stemIt->second.to_chars(), (std::pair<char, char>{'C', 'U'}));
    EXPECT_FLOAT_EQ(stemIt->first, log2f(17328 * (1*600+1) * 2. / 600 / 5 / 103));
    EXPECT_EQ(++stemIt, stem.prio[4].cend());

    // check stem.prio[3] values: {0,0,0,1,0,1,1,0,0,0,0,0,2,0,0,0}
    EXPECT_EQ(stem.prio[3].size(), 4ul);
    EXPECT_FLOAT_EQ(stem.prio[3].cbegin()->first, log2f(17328 * (1*600+1) * 2. / 600 / 5 / 9913));

    EXPECT_EQ(stem.gaps.size(), 5u);
    for (auto const & map : stem.gaps)
    {
        EXPECT_TRUE(map.empty());
    }

    // check the loop
    EXPECT_TRUE(std::holds_alternative<mars::LoopElement>(motif.elements[0]));
    mars::LoopElement loop = std::get<mars::LoopElement>(motif.elements[0]);
    EXPECT_FALSE(loop.is_5prime);
    EXPECT_EQ(loop.prio.size(), 11u);

    // check loop.prio[1] values: {0,0,0,1}
    EXPECT_EQ(loop.prio[1].size(), 1ul);
    auto loopIt = loop.prio[1].cbegin();
    EXPECT_EQ(loopIt->second.to_char(), 'U');
    EXPECT_FLOAT_EQ(loopIt->first, log2f((1*600+1) / 0.3 / 5 / 600));
    EXPECT_EQ(++loopIt, loop.prio[1].cend());

    // check loop.prio[0] values: {0,2,0,3}
    EXPECT_EQ(loop.prio[0].size(), 2ul);
    loopIt = loop.prio[0].cbegin();
    EXPECT_EQ(loopIt->second.to_char(), 'U');
    EXPECT_FLOAT_EQ(loopIt->first, log2f((3*600+1) / 0.3 / 5 / 600));
    ++loopIt;
    EXPECT_EQ(loopIt->second.to_char(), 'C');
    EXPECT_FLOAT_EQ(loopIt->first, log2f((2*600+1) / 0.2 / 5 / 600));
    EXPECT_EQ(++loopIt, loop.prio[0].cend());

    EXPECT_EQ(loop.gaps.size(), 11u);
    for (int idx : {0, 2, 3, 4, 5, 6, 7, 8, 9, 10})
    {
        EXPECT_TRUE(loop.gaps[idx].empty());
    }
    EXPECT_EQ(loop.gaps[1].size(), 1u);
    auto const & gap_entry = loop.gaps[1].cbegin();
    EXPECT_EQ(gap_entry->first, 4);
    EXPECT_EQ(gap_entry->second, 4u);
}
