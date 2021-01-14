#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>

#include <unit/alphabet/semi_alphabet_test_template.hpp>
#include <unit/alphabet/semi_alphabet_constexpr_test_template.hpp>

#include "bi_alphabet.hpp"

TEST(BiAlphabet, Concept)
{
    EXPECT_TRUE(seqan3::writable_semialphabet<mars::bi_alphabet<seqan3::rna4>>);
    EXPECT_TRUE(seqan3::writable_semialphabet<mars::bi_alphabet<seqan3::dna5>>);
    EXPECT_TRUE(seqan3::writable_semialphabet<mars::bi_alphabet<mars::bi_alphabet<seqan3::dna5>>>);
}

TEST(BiAlphabet, AlphabetSize)
{
    EXPECT_EQ(seqan3::alphabet_size<mars::bi_alphabet<seqan3::rna4>>, 16u);
    EXPECT_EQ(seqan3::alphabet_size<mars::bi_alphabet<seqan3::rna5>>, 25u);
}

TEST(BiAlphabet, Construction)
{
    using seqan3::operator""_rna4;

    mars::bi_alphabet<seqan3::rna4> chr1{'A'_rna4, 'C'_rna4};
    // type deduction
    mars::bi_alphabet chr2{'A'_rna4, 'C'_rna4};
    EXPECT_TRUE((std::is_same_v<decltype(chr1), decltype(chr2)>));
    EXPECT_TRUE(chr1 == chr2);
}

TEST(BiAlphabet, Assignment)
{
    using seqan3::operator""_rna5;
    using seqan3::get;

    // default construction
    mars::bi_alphabet<seqan3::rna5> chr{};
    EXPECT_TRUE(chr == (mars::bi_alphabet{'A'_rna5, 'A'_rna5}));

    // assignment
    chr = mars::bi_alphabet{'U'_rna5, 'U'_rna5};
    EXPECT_TRUE(chr == (mars::bi_alphabet{'U'_rna5, 'U'_rna5}));

    // partial assignment
    get<1>(chr) = 'N'_rna5;
    EXPECT_TRUE(chr == (mars::bi_alphabet{'U'_rna5, 'N'_rna5}));
}

TEST(BiAlphabet, GetValue)
{
    using seqan3::operator""_rna5;
    using seqan3::get;

    mars::bi_alphabet chr{'G'_rna5, 'C'_rna5};
    EXPECT_TRUE(get<0>(chr) == 'G'_rna5);
    EXPECT_TRUE(get<1>(chr) == 'C'_rna5);
}

TEST(BiAlphabet, Rank)
{
    using seqan3::operator""_rna4;
    using seqan3::get;

    // retrieve rank
    mars::bi_alphabet chr{'G'_rna4, 'C'_rna4};
    auto rnk = chr.to_rank();
    EXPECT_TRUE((std::is_same_v<decltype(rnk), uint8_t>));
    EXPECT_EQ(rnk, 9);

    // assign rank
    chr.assign_rank(7);
    EXPECT_TRUE(get<0>(chr) == 'C'_rna4);
    EXPECT_TRUE(get<1>(chr) == 'U'_rna4);
}

TEST(BiAlphabet, ToChars)
{
    using seqan3::operator ""_rna4;

    mars::bi_alphabet bi{'G'_rna4, 'C'_rna4};
    auto chrs = bi.to_chars();
    EXPECT_TRUE((std::is_same_v<decltype(chrs), std::pair<char, char>>));
    EXPECT_EQ(chrs.first, 'G');
    EXPECT_EQ(chrs.second, 'C');
}

TEST(BiAlphabet, AssignChars)
{
    using seqan3::operator ""_rna5;
    using seqan3::get;

    mars::bi_alphabet<seqan3::rna5> bi;
    bi.assign_chars('C', 'U');
    EXPECT_TRUE(get<0>(bi) == 'C'_rna5);
    EXPECT_TRUE(get<1>(bi) == 'U'_rna5);
}

TEST(BiAlphabet, CharIsValid)
{
    using seqan3::operator""_rna5;

    // function call from member
    mars::bi_alphabet chr{'G'_rna5, 'C'_rna5};
    EXPECT_TRUE(chr.char_is_valid('U'));
    EXPECT_TRUE(chr.char_is_valid('N'));
    EXPECT_FALSE(chr.char_is_valid('M'));

    // function call from class
    EXPECT_TRUE(mars::bi_alphabet<seqan3::rna4>::char_is_valid('c'));
    EXPECT_FALSE(mars::bi_alphabet<seqan3::rna4>::char_is_valid('N'));
    EXPECT_FALSE(mars::bi_alphabet<seqan3::rna4>::char_is_valid('S'));
}

// Test templates for semi alphabet
using bi_types = ::testing::Types<mars::bi_alphabet<seqan3::rna4>,
                                  mars::bi_alphabet<seqan3::rna5>,
                                  mars::bi_alphabet<seqan3::dna4>,
                                  mars::bi_alphabet<seqan3::dna5>>;

INSTANTIATE_TYPED_TEST_SUITE_P(BiAlphabet, semi_alphabet_test, bi_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(BiAlphabet, semi_alphabet_constexpr, bi_types, );
