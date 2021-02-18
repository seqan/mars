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
    mars::BiDirectionalSearch bds(4);
    EXPECT_NO_THROW(bds.create_index(data("genome.fa")));

    // from archive
    EXPECT_NO_THROW(bds.create_index(data("genome2.fa")));
}

TEST(Index, BiDirectionalSearch)
{
    using seqan3::operator""_rna4;

    mars::BiDirectionalSearch bds(4);
    bds.create_index(data("RF0005.fa"));
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
