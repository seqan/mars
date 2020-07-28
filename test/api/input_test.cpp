#include <gtest/gtest.h>

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <string_view>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "input.hpp"

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(input, read_clustal_file)
{
    std::vector<std::string> names
    {
        {"M83762.1-1031_1093"},
        {"AC008670.6-83725_83795"},
        {"Z82044.1-16031_16103"},
        {"AE004843.1-4972_4900"},
        {"AB042432.1-14140_14072"}
    };

    std::vector<std::vector<seqan3::gapped<seqan3::rna5>>> alignment{5};
    using std::ranges::copy;
    copy(std::string_view{"gcuuuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguagaa-aacu--a-ua---cguuuuaaagcu"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna5>>, std::cpp20::back_inserter(alignment[0]));
    copy(std::string_view{"acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaaau-uuuggugcaacuccaaauaaaagua"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna5>>, std::cpp20::back_inserter(alignment[1]));
    copy(std::string_view{"gcgguuguggcgaag-ugguuaacgcaccagauuguggcucuggcacuc----guggguucgauucccaucaaucgcc"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna5>>, std::cpp20::back_inserter(alignment[2]));
    copy(std::string_view{"gcucauguagc-ucaguugguagagcacacccu----ugguaagggugaggucagcgguucaaauccgcucaugagcu"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna5>>, std::cpp20::back_inserter(alignment[3]));
    copy(std::string_view{"guuucuguagu-ugaau---uacaacgaugauu----uuucaugucauuggu-cgcaguugaaugcuguguagaaaua"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna5>>, std::cpp20::back_inserter(alignment[4]));

    mars::multiple_alignment<seqan3::rna5> msa = mars::read_clustal_file<seqan3::rna5>(data("tRNA.aln"));

    EXPECT_RANGE_EQ(msa.sequences, alignment);
    EXPECT_RANGE_EQ(msa.names, names);
}
