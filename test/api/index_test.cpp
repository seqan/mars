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
    mars::BiDirectionalIndex bds(4);
    EXPECT_NO_THROW(bds.create(data("genome.fa")));
#ifdef SEQAN3_HAS_ZLIB
    std::filesystem::path const indexfile = data("genome.fa.marsindex.gz");
#else
    std::filesystem::path const indexfile = data("genome.fa.marsindex");
#endif
    EXPECT_TRUE(std::filesystem::exists(indexfile));
    std::filesystem::remove(indexfile);

    // from archive
    EXPECT_NO_THROW(bds.create(data("genome2.fa")));

    // from compressed archive
#ifdef SEQAN3_HAS_ZLIB
    EXPECT_NO_THROW(bds.create(data("genome3.fa")));
#endif
}

TEST(Index, BiDirectionalIndex)
{
    using seqan3::operator""_rna4;

    mars::BiDirectionalIndex bds(4);
    bds.create(data("RF0005.fa"));
    mars::bi_alphabet ba{'U'_rna4, 'C'_rna4};

    EXPECT_TRUE(bds.append_loop({1.f, 'A'_rna4}, false));
    EXPECT_TRUE(bds.append_stem({2.f, ba}));
    EXPECT_TRUE(bds.append_loop({0.f, 'A'_rna4}, true));
    EXPECT_TRUE(bds.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_TRUE(bds.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(bds.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(bds.append_stem({0.5f, ba}));
    bds.backtrack();
    EXPECT_TRUE(bds.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(bds.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_FALSE(bds.xdrop());
    EXPECT_TRUE(bds.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_FALSE(bds.append_loop({0.f, 'G'_rna4}, false));

    std::vector<mars::Hit> hits{};
    bds.compute_hits(hits);
    EXPECT_EQ(hits.size(), 5ul);
}
