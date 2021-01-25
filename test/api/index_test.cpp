#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "index.hpp"

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(Index, Create)
{
    // from fasta file
    mars::index_type index1 = mars::create_index(data("genome.fa"));

    // from archive
    mars::index_type index2 = mars::create_index(data("genome2.fa"));

    EXPECT_TRUE(index1 == index2);
}

TEST(Index, BiDirectionalSearch)
{
    using seqan3::operator""_dna4;

    mars::index_type index = mars::create_index(data("RF0005.fa"));
    mars::bi_directional_search bds(index);

    bds.append_3prime('A'_dna4);
    bds.append_stem('U'_dna4, 'C'_dna4);
    bds.append_5prime('A'_dna4);
    bds.append_3prime('G'_dna4);
    bds.append_3prime('A'_dna4);
    bds.append_3prime('A'_dna4);
    bds.append_stem('U'_dna4, 'C'_dna4);
    bds.backtrack();
    bds.append_3prime('A'_dna4);
    bds.append_3prime('G'_dna4);
    bds.append_3prime('G'_dna4);
    bds.append_3prime('G'_dna4);

    size_t num_matches = bds.compute_matches();
    EXPECT_EQ(num_matches, 0ul);

    bds.backtrack();
    num_matches = bds.compute_matches();
    EXPECT_EQ(num_matches, 5ul);
}
