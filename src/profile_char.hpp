#pragma once

#include <array>
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/utility/views/zip.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/array.hpp>
#endif

namespace mars
{

/*!
 * \brief Stores the expected background distribution for structured RNA sequences.
 */
struct
{
    std::array<float, 16> const r16
        {
            // Source:
            // Olson, W. K., Esguerra, M., Xin, Y., & Lu, X. J. (2009).
            // New information content in RNA base pairing deduced from quantitative analysis of high-resolution structures.
            // Methods (San Diego, Calif.), 47(3), 177–186. https://doi.org/10.1016/j.ymeth.2008.12.003

            log2f( 384.f     / 17328), // AA
            log2f( 313.f / 2 / 17328), // AC
            log2f( 980.f / 2 / 17328), // AG
            log2f(3975.f / 2 / 17328), // AU
            log2f( 313.f / 2 / 17328), // CA
            log2f(  63.f     / 17328), // CC
            log2f(9913.f / 2 / 17328), // CG
            log2f( 103.f / 2 / 17328), // CU
            log2f( 980.f / 2 / 17328), // GA
            log2f(9913.f / 2 / 17328), // GC
            log2f( 128.f     / 17328), // GG
            log2f(1282.f / 2 / 17328), // GU
            log2f(3975.f / 2 / 17328), // UA
            log2f( 103.f / 2 / 17328), // UC
            log2f(1282.f / 2 / 17328), // UG
            log2f( 187.f     / 17328)  // UU
        };

    std::array<float, 4> const r4
        {
            log2f(0.3f), // A
            log2f(0.2f), // C
            log2f(0.2f), // G
            log2f(0.3f)  // U
        };

    template <unsigned short size> requires (size == 4 || size == 16)
    std::array<float, size> const & get() const
    {
        if constexpr (size == 4)
            return r4;
        else
            return r16;
    }
} BackgroundDistribution;

/*!\interface BiAlphabetConcept <>
 * \brief A concept that checks whether t is a Bialphabet.
 * \tparam t The type to be checked.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT BiAlphabetConcept =
requires (t v)
{
    { seqan3::alphabet_size<t> };
    { v.assign_chars('c', 'c') };
    { v.to_chars() };
};
//!\endcond

/*!
 * \brief Stores the frequency of characters at a specific position.
 * \tparam alph_type The legal alphabet for the profile entries.
 */
template <seqan3::semialphabet alph_type>
class profile_char
{
private:
    //! \brief The number of profile entries, i.e. the size of the underlying alphabet.
    static constexpr size_t const size{seqan3::alphabet_size<alph_type>};

    //! \brief The internal representation of a single count.
    static constexpr uint32_t const one{600};

    //! \brief Convert wildcard characters into their components.
    static std::string compose(char chr)
    {
        switch (chr)
        {
            case 'A': return "A";
            case 'C': return "C";
            case 'G': return "G";
            case 'U': return "U";
            case 'T': return "T";
            case 'M': return "AC";
            case 'R': return "AG";
            case 'W': return "AU";
            case 'Y': return "CU";
            case 'S': return "CG";
            case 'K': return "GU";
            case 'V': return "ACG";
            case 'H': return "ACU";
            case 'D': return "AGU";
            case 'B': return "CGU";
            case 'N':
                if constexpr (!alph_type::char_is_valid('N')) // split 'N' if not supported
                    return "ACGU";
                else
                    return "N";
            default: throw std::runtime_error(std::string{"Invalid character found: "} + chr);
        }
    }

    //! \brief The number of occurrences of each character.
    std::array<uint32_t, size> tally{};

public:
    /*!
     * \name Constructors, destructor and assignment
     * \{
     */
    constexpr profile_char()                                 noexcept = default; //!< Defaulted.
    constexpr profile_char(profile_char const &)             noexcept = default; //!< Defaulted.
    constexpr profile_char(profile_char &&)                  noexcept = default; //!< Defaulted.
    constexpr profile_char & operator=(profile_char const &) noexcept = default; //!< Defaulted.
    constexpr profile_char & operator=(profile_char &&)      noexcept = default; //!< Defaulted.
    ~profile_char()                                          noexcept = default; //!< Defaulted.

    /*!
     * \brief Increase the character count by 1.
     * \param rnk The rank of the character of which the count is incremented.
     */
    void increment(seqan3::alphabet_rank_t<alph_type> rnk)
    {
        assert(rnk < size);
        tally[rnk] += one;
    }

    /*!
     * \brief Increase the character count by 1.
     * \param chr The character of which the count is incremented.
     */
    void increment(alph_type chr)
    {
        tally[chr.to_rank()] += one;
    }

    /*!
     * \brief Increase the character count by 1.
     * \tparam ext_alph_type A compatible nucleotide alphabet type that may be smaller than `alph_type`.
     * \param chr The character of which the count is incremented.
     */
    template <seqan3::nucleotide_alphabet ext_alph_type>
    //!\cond
        requires seqan3::nucleotide_alphabet<alph_type> &&
                 (seqan3::alphabet_size<ext_alph_type> <= seqan3::alphabet_size<alph_type>)
    //!\endcond
    void increment(ext_alph_type chr)
    {
        increment(alph_type{chr}); // DNA-RNA conversion or convert into a larger alphabet
    }

    /*!
     * \brief Increase the character count (by 1 in total).
     * \tparam ext_alph_type The extended alphabet type that may contain wildcards; must be a nucleotide alphabet.
     * \param chr The character of which the count is incremented.
     *
     * \details
     * If a wildcard is given, the counts of all matching characters are increased by the same fraction
     * (e.g. M = 1/2 A + 1/2 C).
     */
    template <seqan3::nucleotide_alphabet ext_alph_type>
    //!\cond
        requires seqan3::writable_alphabet<alph_type> &&
                 (seqan3::alphabet_size<ext_alph_type> > seqan3::alphabet_size<alph_type>)
    //!\endcond
    void increment(ext_alph_type chr)
    {
        std::string composition = compose(chr.to_char());
        for (char x : composition)
            tally[alph_type{}.assign_char(x).to_rank()] += one/composition.size();
    }

    /*!
     * \brief This is an overload for gapped alphabets. Increase the character count by 1 unless it is a gap.
     * \tparam inner_alph_type The underlying alphabet type without the gap; must be a writable alphabet.
     * \param chr The character of which the count is incremented.
     * \returns True, if `chr` is a gap character, false otherwise.
     */
    template <seqan3::writable_alphabet inner_alph_type>
    bool increment(seqan3::gapped<inner_alph_type> chr)
    {
        if (chr == seqan3::gap())
            return true;

        increment(inner_alph_type{}.assign_char(chr.to_char()));
        return false;
    }

    /*!
     * \brief Increase the character count (by 1 in total) for bi-alphabets.
     * \tparam ext_alph_type The extended alphabet type that may contain wildcards; must be a nucleotide alphabet.
     * \param chr1 The first character of the bi-character.
     * \param chr2 The second character of the bi-character.
     *
     * \details
     * If a wildcard is given, the counts of all matching characters are increased by the same fraction
     * (e.g. M = 1/2 A + 1/2 C).
     */
    template <seqan3::nucleotide_alphabet ext_alph_type>
    //!\cond
        requires BiAlphabetConcept<alph_type>
    //!\endcond
    void increment(ext_alph_type chr1, ext_alph_type chr2)
    {
        std::string composition1 = compose(seqan3::to_char(chr1));
        std::string composition2 = compose(seqan3::to_char(chr2));
        size_t len = composition1.size() * composition2.size();
        for (char x1 : composition1)
            for (char x2 : composition2)
                tally[alph_type{}.assign_chars(x1, x2).to_rank()] += one/len;
    }

    /*!
     * \brief This is an overload for bi-alphabets with gaps. Increase the character count by 1 unless it has a gap.
     * \tparam inner_alph_type The underlying alphabet type without the gap; must be a writable alphabet.
     * \param chr1 The first (gapped) character of the bi-character.
     * \param chr2 The second (gapped) character of the bi-character.
     * \returns True, if `chr1` or `chr2` is a gap character, false otherwise.
     */
    template <seqan3::writable_alphabet inner_alph_type>
    //!\cond
        requires BiAlphabetConcept<alph_type>
    //!\endcond
    bool increment(seqan3::gapped<inner_alph_type> chr1, seqan3::gapped<inner_alph_type> chr2)
    {
        if (chr1 == seqan3::gap() || chr2 == seqan3::gap())
            return true;

        increment(inner_alph_type{}.assign_char(chr1.to_char()), inner_alph_type{}.assign_char(chr2.to_char()));
        return false;
    }

    /*!
     * \brief Retrieve the quantity of a character.
     * \param chr The character of interest.
     * \return The quantity of the character in the profile.
     *
     * \note
     *  The returned number can be non-integral if wildcards were present.
     */
    float quantity(alph_type chr) const
    {
        return 1.f * tally[chr.to_rank()] / one;
    }

    /*!
     * \brief Retrieve the quantity of a character.
     * \param rank The rank of the character of interest.
     * \return The quantity of the character in the profile.
     *
     * \note
     *  The returned number can be non-integral if wildcards were present.
     */
    [[nodiscard]] float quantity(uint32_t rank) const
    {
        assert(rank < size);
        return tally[rank] / static_cast<float>(one);
    }

    /*!
     * \brief Retrieve the whole character profile.
     * \return An array that contains the quantities of the profile characters in alphabetical order.
     *
     * \note
     *  The contained numbers can be non-integral if wildcards were present.
     */
    std::array<float, size> quantities() const
    {
        std::array<float, size> tmp;
        std::transform(tally.begin(), tally.end(), tmp.begin(), [] (uint32_t x)
        {
            return x / static_cast<float>(one);
        });
        return tmp;
    }

    /*!
     * \brief Retrieve the log quantities relative to the background distribution.
     * \param depth The number of sequences used to create the profile (for normalization).
     * \return A priority queue with logarithmic scores and the respective RNA characters.
     */
    std::set<std::pair<float, alph_type>> priority(size_t depth) const
    {
        std::set<std::pair<float, alph_type>> result{};
        for (auto && [idx, bg] : seqan3::views::zip(std::ranges::views::iota(0), BackgroundDistribution.get<size>()))
            result.emplace(log2f((tally[idx] + 1.) / one / depth) - bg, alph_type{}.assign_rank(idx));
        return result;
    }

#if SEQAN3_WITH_CEREAL
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(tally);
    }
#endif
};

/*!
 * \brief Stream a representation of a character profile.
 * \tparam alph_type The alphabet type of the profile.
 * \tparam ostream_type The stream type.
 * \param os The stream where the representation is appended.
 * \param chr The character profile that should be printed.
 * \return The output stream.
 * \relates profile_char
 */
template <seqan3::semialphabet alph_type, typename ostream_type>
inline ostream_type & operator<<(ostream_type & os, profile_char<alph_type> const & chr)
{
    os << "(" << chr.quantity(0);
    for (size_t idx = 1; idx < seqan3::alphabet_size<alph_type>; ++idx)
    {
        os << "," << chr.quantity(idx);
    }
    os << ")";
    return os;
}

} // namespace mars
