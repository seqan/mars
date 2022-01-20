#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/rna4.hpp>

#include "bi_alphabet.hpp"
#include "index.hpp"
#include "search.hpp"
#include "settings.hpp"

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(Index, Create)
{
    // from fasta file
    mars::BiDirectionalIndex bds{};
    mars::settings.genome_file = data("genome.fa");
    mars::settings.compress_index = true;
    EXPECT_NO_THROW(bds.create());
#ifdef SEQAN3_HAS_ZLIB
    std::filesystem::path const indexfile = data("genome.fa.marsindex.gz");
#else
    std::filesystem::path const indexfile = data("genome.fa.marsindex");
#endif
    EXPECT_TRUE(std::filesystem::exists(indexfile));
    std::filesystem::remove(indexfile);

    // from archive
    mars::settings.genome_file = data("genome2.fa");
    EXPECT_NO_THROW(bds.create());

    // from compressed archive
#ifdef SEQAN3_HAS_ZLIB
    mars::settings.genome_file = data("genome3.fa");
    EXPECT_NO_THROW(bds.create());
#endif
}

TEST(Index, BiDirectionalIndex)
{
    using seqan3::operator""_rna4;

    mars::BiDirectionalIndex bds{};
    mars::settings.genome_file = data("RF00005.fa");
    bds.create();
    mars::bi_alphabet ba{'U'_rna4, 'C'_rna4};
    mars::HitStore hits(10);
    mars::StemloopMotif motif{0, {27, 47}};
    mars::SearchInfo info{bds.raw(), motif, hits};

    EXPECT_TRUE(info.append_loop({1.f, 'A'_rna4}, false));
    EXPECT_TRUE(info.append_stem({2.f, ba}));
    EXPECT_TRUE(info.append_loop({0.f, 'A'_rna4}, true));
    EXPECT_TRUE(info.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_TRUE(info.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(info.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(info.append_stem({0.5f, ba}));
    info.backtrack();
    EXPECT_TRUE(info.append_loop({0.f, 'A'_rna4}, false));
    EXPECT_TRUE(info.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_FALSE(info.xdrop());
    EXPECT_TRUE(info.append_loop({0.f, 'G'_rna4}, false));
    EXPECT_FALSE(info.append_loop({0.f, 'G'_rna4}, false));

    EXPECT_NO_THROW(info.compute_hits());
}
