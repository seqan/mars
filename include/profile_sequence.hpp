#pragma once

#include <array>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

namespace mars
{

/*!\brief Stores the frequency of characters at a specific position.
 * \tparam alph_type The legal alphabet for the profile entries.
 */
template <seqan3::alphabet alph_type>
class profile_char
{
private:
    //!\brief The number of profile entries, i.e. the size of the underlying alphabet.
    static constexpr size_t size{seqan3::alphabet_size<alph_type>};

    //!\brief The internal representation of a single count.
    static constexpr uint32_t one{600};

    //!\brief We store the count values here.
    std::array<uint32_t, size> tally;

    /*!\brief Helper function to increase the respective character frequency by a fraction.
     * \tparam div The fraction of a count, e.g. 3 increases the frequency by 1/3.
     * \param x The character of which the count is increased.
     */
    template <uint32_t div>
    void incr (char x)
    {
        tally[alph_type{}.assign_char(x).to_rank()] += one/div;
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr profile_char()                                 noexcept = default; //!< Defaulted.
    constexpr profile_char(profile_char const &)             noexcept = default; //!< Defaulted.
    constexpr profile_char(profile_char &&)                  noexcept = default; //!< Defaulted.
    constexpr profile_char & operator=(profile_char const &) noexcept = default; //!< Defaulted.
    constexpr profile_char & operator=(profile_char &&)      noexcept = default; //!< Defaulted.
    ~profile_char()                                          noexcept = default; //!< Defaulted.

    /*!\brief Increase the character count by 1.
     * \param chr The character of which the count is incremented.
     */
    void increment(alph_type chr)
    {
        tally[chr.to_rank()] += one;
    }

    /*!\brief Increase the character count (by 1 in total).
     * \tparam ext_alph_type The extended alphabet type that may contain wildcards; must be a nucleotide alphabet.
     * \param chr The character of which the count is incremented.
     *
     * \details
     * If a wildcard is given, the counts of all matching characters are increased by the same fraction
     * (e.g. M = 1/2 A + 1/2 C).
     */
    template <seqan3::nucleotide_alphabet ext_alph_type>
    void increment(ext_alph_type chr)
    {
        static_assert(seqan3::nucleotide_alphabet<alph_type>,
                      "Wildcard conversion is only implemented for nucleotide alphabets.");

        if constexpr (seqan3::alphabet_size<ext_alph_type> <= seqan3::alphabet_size<alph_type>)
        {
            increment(alph_type{chr}); // DNA-RNA conversion or convert into a larger alphabet
        }
        else
        {
            switch (chr.to_char())
            {
                case 'M': incr<2>('A'); incr<2>('C'); break;
                case 'R': incr<2>('A'); incr<2>('G'); break;
                case 'W': incr<2>('A'); incr<2>('T'); break;
                case 'Y': incr<2>('C'); incr<2>('T'); break;
                case 'S': incr<2>('C'); incr<2>('G'); break;
                case 'K': incr<2>('G'); incr<2>('T'); break;
                case 'V': incr<3>('A'); incr<3>('C'); incr<3>('G'); break;
                case 'H': incr<3>('A'); incr<3>('C'); incr<3>('T'); break;
                case 'D': incr<3>('A'); incr<3>('G'); incr<3>('T'); break;
                case 'B': incr<3>('C'); incr<3>('G'); incr<3>('T'); break;
                case 'N':
                    if constexpr (!seqan3::char_is_valid_for<alph_type>('N')) // split 'N' if not supported
                    {
                        incr<4>('A'); incr<4>('C'); incr<4>('G'); incr<4>('T');
                        break;
                    }
                default: increment(alph_type{chr});
            }
        }
    }

    /*!\brief Retrieve the quantity of a character.
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

    /*!\brief Retrieve the quantity of a character.
     * \param chr The rank of the character of interest.
     * \return The quantity of the character in the profile.
     *
     * \note
     *  The returned number can be non-integral if wildcards were present.
     */
    [[nodiscard]] float quantity(uint32_t rank) const
    {
        assert(rank < size);
        return 1.f * tally[rank] / one;
    }

    /*!\brief Retrieve the whole character profile.
     * \return An array that contains the quantities of the profile characters in alphabetical order.
     *
     * \note
     *  The contained numbers can be non-integral if wildcards were present.
     */
    std::array<float, size> quantities() const
    {
        std::array<float, size> tmp;
        std::transform(tally.begin(), tally.end(), tmp.begin(), [] (uint32_t x) { return 1.f * x / one; });
        return std::move(tmp);
    }
};

/*!\brief Stream a representation of a character profile.
 * \tparam alph_type The alphabet type of the profile.
 * \tparam ostream_type The stream type.
 * \param os The stream where the representation is appended.
 * \param chr The character profile that should be printed.
 * \return The output stream.
 */
template <seqan3::alphabet alph_type, typename ostream_type>
inline ostream_type & operator<<(ostream_type & os, profile_char<alph_type> const & chr)
{
    os << "(" << alph_type{}.assign_rank(0).to_char() << ":" << chr.quantity(0);
    for (size_t idx = 1; idx < seqan3::alphabet_size<alph_type>; ++idx)
    {
        os << "," << alph_type{}.assign_rank(idx).to_char() << ":" << chr.quantity(idx);
    }
    os << ")";
    return os;
}

} // namespace mars
