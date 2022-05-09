#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "bi_alphabet.hpp"
#include "profile_char.hpp"

TEST(ProfileChar, SimpleIncrementAndQuantity)
{
    using seqan3::operator""_rna4;
    mars::profile_char<seqan3::rna4> prof{};
    prof.increment('U'_rna4);
    prof.increment('A'_rna4);
    prof.increment(0); // A (Rna4)
    EXPECT_FLOAT_EQ(prof.quantity('A'_rna4), 2);
    EXPECT_FLOAT_EQ(prof.quantity('C'_rna4), 0);
    EXPECT_FLOAT_EQ(prof.quantity('G'_rna4), 0);
    EXPECT_FLOAT_EQ(prof.quantity('U'_rna4), 1);
    EXPECT_FLOAT_EQ(prof.quantity(0), 2);
    EXPECT_FLOAT_EQ(prof.quantity(1), 0);
    EXPECT_FLOAT_EQ(prof.quantity(2), 0);
    EXPECT_FLOAT_EQ(prof.quantity(3), 1); //                 A  C  G  U
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 4>{1200,0,0,600}));
}

TEST(ProfileChar, ConvertIncrementRna15Rna4)
{
    using seqan3::operator""_rna15;
    mars::profile_char<seqan3::rna4> prof{};
    prof.increment('T'_rna15); // U
    prof.increment('A'_rna15);
    prof.increment('N'_rna15); // ACGU
    prof.increment('N'_rna15); // ACGU
    prof.increment('M'_rna15); // AC
    prof.increment('S'_rna15); // CG                         A  C    G  U
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 4>{1200,900,600,900}));

    prof.increment('V'_rna15); // ACG
    prof.increment('H'_rna15); // ACU
    prof.increment('D'_rna15); // AGU
    prof.increment('B'_rna15); // CGU                        A  C    G  U
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 4>{1800,1500,1200,1500}));
}

TEST(ProfileChar, ConvertIncrementDna15Dna5)
{
    using seqan3::operator""_dna15;
    mars::profile_char<seqan3::dna5> prof{};
    prof.increment('U'_dna15); // T
    prof.increment('A'_dna15);
    prof.increment('N'_dna15); // N
    prof.increment('N'_dna15); // N
    prof.increment('M'_dna15); // AC
    prof.increment('S'_dna15); // CG                         A    C  G    N  U
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 5>{900,600,300,1200,600}));

    prof.increment('R'_dna15); // AG
    prof.increment('W'_dna15); // AU
    prof.increment('Y'_dna15); // CU
    prof.increment('K'_dna15); // GU                         A    C    G    N  U
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 5>{1500,900,900,1200,1500}));
}

TEST(ProfileChar, ConvertIncrementRna4Rna15)
{
    using seqan3::operator""_rna4;
    mars::profile_char<seqan3::rna15> prof{};
    prof.increment('T'_rna4); // U
    prof.increment('A'_rna4);
    prof.increment(4);        // G (Rna15)
    prof.increment('U'_rna4);
    prof.increment('G'_rna4);
    prof.increment('C'_rna4); //                              A  B  C  D  G  H  K  M  N  R  S  U  V  W  Y
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 15>{600,0,600,0,1200,0,0,0,0,0,0,1200,0,0,0}));
}

TEST(ProfileChar, GappedAlphabet)
{
    using seqan3::operator""_rna15;
    mars::profile_char<seqan3::rna4> prof{};

    seqan3::gapped<seqan3::rna15> chr{'m'_rna15};
    EXPECT_FALSE(prof.increment(chr));

    chr = seqan3::gap();
    EXPECT_TRUE(prof.increment(chr));
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 4>{300,300,0,0}));
}

TEST(ProfileChar, StreamOperator)
{
    using seqan3::operator""_rna4;
    mars::profile_char<seqan3::rna4> prof{};
    prof.increment('U'_rna4);
    prof.increment('A'_rna4);
    prof.increment('A'_rna4);

    std::ostringstream os{};
    os << prof;
    EXPECT_STREQ(os.str().c_str(), "(2,0,0,1)");
}

TEST(ProfileChar, BiAlphabet)
{
    using seqan3::operator""_rna4;
    using seqan3::operator""_rna15;

    mars::profile_char<mars::bi_alphabet<seqan3::rna4>> prof{};
    prof.increment(2);
    prof.increment(4);
    prof.increment(6); //                                     AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 16>{0,0,600,0,600,0,600,0,0,0,0,0,0,0,0,0}));

    // Assign with alphabet_type.
    prof.increment({'C'_rna4, 'G'_rna4});
    prof.increment({'U'_rna4, 'A'_rna4});
    prof.increment({'U'_rna4, 'U'_rna4});
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 16>{0,0,600,0,600,0,1200,0,0,0,0,0,600,0,0,600}));

    // Assign with two arguments.
    prof.increment('C'_rna4, 'G'_rna4);
    prof.increment('A'_rna4, 'A'_rna4);
    prof.increment('U'_rna4, 'U'_rna4);
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 16>{600,0,600,0,600,0,1800,0,0,0,0,0,600,0,0,1200}));

    // only N
    mars::profile_char<mars::bi_alphabet<seqan3::rna5>> n5{};
    n5.increment('N'_rna15, 'N'_rna15);
    EXPECT_RANGE_EQ(n5.quantities(), (std::array<float, 25>{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,600,0,0,0,0,0,0}));
    n5.increment('C'_rna4, 'C'_rna4);
    EXPECT_RANGE_EQ(n5.quantities(), (std::array<float, 25>{0,0,0,0,0,0,600,0,0,0,0,0,0,0,0,0,0,0,600,0,0,0,0,0,0}));
}

TEST(ProfileChar, BiAlphabetGaps)
{
    using seqan3::operator""_rna15;
    using seqan3::operator""_rna5;

    seqan3::gapped<seqan3::rna15> r{'R'_rna15};
    seqan3::gapped<seqan3::rna15> n{'N'_rna15};
    seqan3::gapped<seqan3::rna15> a{'A'_rna15};
    seqan3::gapped<seqan3::rna15> g{seqan3::gap()};

    mars::profile_char<mars::bi_alphabet<seqan3::gapped<seqan3::rna4>>> prof{};
    prof.increment(r, r);
    prof.increment(g, n);
    prof.increment(g, g);
    prof.increment(r, r);
    prof.increment(a, a); // AA AC AG AU A- CA CC CG CU C- GA GC GG GU G- UA UC UG UU U- -A -C -G -U --
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 25>{900,0,300,0,0,0,0,0,0,0,300,0,300,0,0,0,0,0,0,0,
                                                              150,150,150,150,600}));

    prof.increment('A'_rna15, 'N'_rna15);
    prof.increment('A'_rna15, 'N'_rna15);
    prof.increment('C'_rna15, 'C'_rna15);
    EXPECT_RANGE_EQ(prof.quantities(), (std::array<float, 25>{1200,300,600,300,0,0,600,0,0,0,300,0,300,0,0,0,0,0,0,0,
                                                              150,150,150,150,600}));
}
