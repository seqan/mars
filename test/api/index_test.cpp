#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/rna4.hpp>

#include "index.hpp"

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(Index, Create)
{
    // from fasta file
    mars::Index index1 = mars::create_index(data("genome.fa"));

    // from archive
    mars::Index index2 = mars::create_index(data("genome2.fa"));

    EXPECT_TRUE(index1 == index2);
}

TEST(Index, BiDirectionalSearch)
{
    using seqan3::operator""_rna4;

    mars::Index index = mars::create_index(data("RF0005.fa"));
    mars::BiDirectionalSearch bds(index, 4);
    mars::bi_alphabet ba{'U'_rna4, 'C'_rna4};

    bds.append_loop({0.f, 'A'_rna4}, false);
    bds.append_stem({0.f, ba});
    bds.append_loop({0.f, 'A'_rna4}, true);
    bds.append_loop({0.f, 'G'_rna4}, false);
    bds.append_loop({0.f, 'A'_rna4}, false);
    bds.append_loop({0.f, 'A'_rna4}, false);
    bds.append_stem({0.f, ba});
    bds.backtrack();
    bds.append_loop({0.f, 'A'_rna4}, false);
    bds.append_loop({0.f, 'G'_rna4}, false);
    bds.append_loop({0.f, 'G'_rna4}, false);
    bds.append_loop({0.f, 'G'_rna4}, false);

    size_t num_matches = bds.compute_matches();
    EXPECT_EQ(num_matches, 0ul);

    bds.backtrack();
    num_matches = bds.compute_matches();
    EXPECT_EQ(num_matches, 5ul);
}
