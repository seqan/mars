#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>

#include <unit/alphabet/semi_alphabet_test_template.hpp>
#include <unit/alphabet/semi_alphabet_constexpr_test_template.hpp>

#include "bi_alphabet.hpp"

TEST(BiAlphabet, AlphabetSize)
{
    using alph_t = mars::bi_alphabet<seqan3::rna4>;
    EXPECT_EQ(seqan3::alphabet_size<alph_t>, 16u);
}

TEST(BiAlphabet, Construction)
{
    using seqan3::operator""_rna4;
    using alph_t = mars::bi_alphabet<seqan3::rna4>;
    alph_t chr{'A'_rna4, 'C'_rna4};
    alph_t chr2 = chr;
    EXPECT_EQ(chr2, chr);
    EXPECT_EQ(chr2, (alph_t{'A'_rna4, 'C'_rna4}));
    alph_t chr3{'G'_rna4};
    EXPECT_EQ(chr3, (alph_t{'G'_rna4, 'G'_rna4}));
}

TEST(BiAlphabet, Assignment)
{
    using seqan3::operator""_rna5;
    using seqan3::get;
    using alph_t = mars::bi_alphabet<seqan3::rna5>;
    alph_t chr{};
    EXPECT_EQ(chr, (alph_t{'A'_rna5, 'A'_rna5}));
    chr = 'U'_rna5;
    EXPECT_EQ(chr, (alph_t{'U'_rna5, 'U'_rna5}));
    get<1>(chr) = 'N'_rna5;
    EXPECT_EQ(chr, (alph_t{'U'_rna5, 'N'_rna5}));
}

TEST(BiAlphabet, Value)
{
    using seqan3::operator""_rna5;
    using seqan3::get;
    mars::bi_alphabet<seqan3::rna5> chr{'G'_rna5, 'C'_rna5};
    EXPECT_EQ(get<0>(chr), 'G'_rna5);
    EXPECT_EQ(get<1>(chr), 'C'_rna5);
}

using bi_types = ::testing::Types<mars::bi_alphabet<seqan3::rna4>,
                                  mars::bi_alphabet<seqan3::rna5>,
                                  mars::bi_alphabet<seqan3::dna4>,
                                  mars::bi_alphabet<seqan3::dna5>>;

INSTANTIATE_TYPED_TEST_SUITE_P(BiAlphabet, semi_alphabet_test, bi_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(BiAlphabet, semi_alphabet_constexpr, bi_types, );
