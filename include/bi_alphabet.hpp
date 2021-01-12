#pragma once

#include <array>
#include <vector>

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace mars
{

//!\brief Only seqan3::get works with bi_alphabet.
using seqan3::get;

/*!\brief A seqan3::alphabet_tuple_base that joins an alphabet with itself to represent a double character.
 * \tparam alphabet_t Type of the underlying alphabet; must model seqan3::writable_semialphabet.
 *
 * \details
 * This seqan3::alphabet_tuple_base itself models seqan3::writable_semialphabet.
 */
template <seqan3::writable_semialphabet alphabet_t>
class bi_alphabet : public seqan3::alphabet_tuple_base<bi_alphabet<alphabet_t>, alphabet_t, alphabet_t>
{
private:
    //!\brief The base type.
    using base_type = seqan3::alphabet_tuple_base<bi_alphabet<alphabet_t>, alphabet_t, alphabet_t>;

public:
    //!\brief The template parameter as member type.
    using alphabet_type = alphabet_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr bi_alphabet()                                 noexcept = default; //!< Defaulted.
    constexpr bi_alphabet(bi_alphabet const &)              noexcept = default; //!< Defaulted.
    constexpr bi_alphabet(bi_alphabet &&)                   noexcept = default; //!< Defaulted.
    constexpr bi_alphabet & operator =(bi_alphabet const &) noexcept = default; //!< Defaulted.
    constexpr bi_alphabet & operator =(bi_alphabet &&)      noexcept = default; //!< Defaulted.
    ~bi_alphabet()                                          noexcept = default; //!< Defaulted.

    using base_type::base_type; // Inherit non-default constructors
    //!\}

    // Inherit operators from base
    using base_type::operator=;

    //!\brief Validate whether a character is valid in the underlying alphabet.
    static constexpr bool char_is_valid(char c) noexcept
    {
        return seqan3::char_is_valid_for<alphabet_type>(c);
    }
};

//!\brief Type deduction guide enables usage of bi_alphabet without specifying template args.
//!\relates bi_alphabet
template <seqan3::writable_semialphabet alphabet_type>
bi_alphabet(alphabet_type &&, alphabet_type &&) -> bi_alphabet<std::decay_t<alphabet_type>>;

} // namespace mars
