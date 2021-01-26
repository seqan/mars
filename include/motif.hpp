#pragma once

#include <ostream>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna4.hpp>

#include "bi_alphabet.hpp"
#include "multiple_alignment.hpp"
#include "profile_char.hpp"

namespace mars
{

//! \brief Store {min, mean, max} of a distribution.
typedef std::tuple<uint16_t, float, uint16_t> stat_type;

//! \brief The boundaries of a stemloop.
typedef std::pair<size_t, size_t> stemloop_type;

//! \brief A loop element in a stemloop.
struct loop_element
{
    stat_type length;
    std::vector<profile_char<seqan3::rna4>> profile;
    std::vector<std::unordered_map<uint16_t, uint16_t>> gaps;
    bool is_5prime;
};

//! \brief A stem element in a stemloop.
struct stem_element
{
    stat_type length;
    std::vector<profile_char<bi_alphabet<seqan3::rna4>>> profile;
    std::vector<std::unordered_map<uint16_t, uint16_t>> gaps;
};

//! \brief A stemloop motif consists of a series of loop and stem elements.
struct stemloop_motif
{
    //! \brief The position interval, where the stemloop is located in the alignment.
    stemloop_type bounds;

    //! \brief The length statistics of the stemloop.
    stat_type length;

    //! \brief A vector of loop and stem elements that the stemloop consists of.
    std::vector<std::variant<loop_element, stem_element>> elements;

    /*!
     * \brief Add a new stem to this motif.
     * \return a reference to the new stem element.
     */
    stem_element & new_stem();

    /*!
     * \brief Add a new loop to this motif.
     * \param is_5prime True for 5' loops, false for 3' loops.
     * \return a reference to the new loop element.
     */
    loop_element & new_loop(bool is_5prime);
};

/*!
 * \brief Stream a representation of a stem loop motif.
 * \param os The output stream.
 * \param motif The motif to be printed.
 * \return the stream with the motif representation appended.
 */
std::ostream & operator<<(std::ostream & os, stemloop_motif const & motif);

/*!
 * \brief Extract the positions of the stem loops.
 * \param bpseq The base pairing at each position.
 * \param plevel The pseudoknot level at each position.
 * \return a vector of stemloop positions.
 */
std::vector<stemloop_type> detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel);

/*!
 * \brief Create a motif from a given stemloop.
 * \param msa The multiple structural alignment.
 * \param bpseq The base pairing at each position.
 * \param pos The stem loop position.
 * \return the motif.
 */
stemloop_motif analyze_stem_loop(msa_type const & msa,
                                  std::vector<int> const & bpseq,
                                  stemloop_type const & pos);

} // namespace mars
