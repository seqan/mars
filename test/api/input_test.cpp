#include <gtest/gtest.h>

#include <seqan3/std/filesystem>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <string_view>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "multiple_alignment.hpp"

// Generate the full path of a test input file that is provided in the data directory.
std::filesystem::path data(std::string const & filename)
{
    return std::filesystem::path{std::string{DATADIR}}.concat(filename);
}

TEST(ClustalInput, ReadFile)
{
    std::vector<std::string> names
    {
        {"M83762.1-1031_1093"},
        {"AC008670.6-83725_83795"},
        {"Z82044.1-16031_16103"},
        {"AE004843.1-4972_4900"},
        {"AB042432.1-14140_14072"}
    };

    std::vector<std::vector<seqan3::gapped<seqan3::rna15>>> alignment{5};
    using std::ranges::copy;
    copy(std::string_view{"gcuuuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguagaa-aacu--a-ua---cguuuuaaagcu"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(alignment[0]));
    copy(std::string_view{"acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaaau-uuuggugcaacuccaaauaaaagua"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(alignment[1]));
    copy(std::string_view{"gcgguuguggcgaag-ugguuaacgcaccagauuguggcucuggcacuc----guggguucgauucccaucaaucgcc"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(alignment[2]));
    copy(std::string_view{"gcucauguagc-ucaguugguagagcacacccu----ugguaagggugaggucagcgguucaaauccgcucaugagcu"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(alignment[3]));
    copy(std::string_view{"guuucuguagu-ugaau---uacaacgaugauu----uuucaugucauuggu-cgcaguugaaugcuguguagaaaua"}
        | seqan3::views::char_to<seqan3::gapped<seqan3::rna15>>, std::cpp20::back_inserter(alignment[4]));

    mars::Msa msa = mars::read_msa(data("tRNA.aln"));

    EXPECT_RANGE_EQ(msa.sequences, alignment);
    EXPECT_RANGE_EQ(msa.names, names);
}

TEST(ClustalInput, FailFileNotFound)
{
    EXPECT_THROW(mars::read_msa(std::filesystem::path{"not_exist.aln"}),
                 seqan3::file_open_error);
}

TEST(ClustalInput, FailClustalHeader)
{
    std::stringstream str{"CLUSTER FORMAT\n\n"};
    EXPECT_THROW(mars::read_msa(str), seqan3::parse_error);
}

TEST(ClustalInput, FailWrongSequenceId)
{
    std::stringstream str{"CLUSTAL FORMAT\n"
                          "\n"
                          "M83762.1-1031_1093      gcuuuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguag\n"
                          "AC008670.6-83725_83795  acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaa\n"
                          "                                 *               *                        \n"
                          "\n"
                          "M83762.1-1031_1093      aa-aacu--a-ua---cguuuuaaagcu\n"
                          "wrong-sequence-id       au-uuuggugcaacuccaaauaaaagua\n"
                          "                                    *               \n"
                          "\n"};
    EXPECT_THROW(mars::read_msa(str), seqan3::parse_error);
}

TEST(ClustalInput, FailInvalidCharacter)
{
    std::stringstream str{"CLUSTAL FORMAT\n"
                          "\n"
                          "M83762.1-1031_1093      gcuzuaaaagc-uuu---gcugaagcaacggcc----uuguaagucguag\n"
                          "AC008670.6-83725_83795  acuuuuaaagg-aua-acagccauccguugguc----uuaggccccaaaa\n"
                          "                                 *               *                        \n"
                          "\n"};
    EXPECT_THROW(mars::read_msa(str), seqan3::parse_error);
}

