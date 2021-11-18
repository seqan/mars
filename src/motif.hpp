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

//! \brief Type for motif id. We expect to find less than 256 motifs.
using MotifNum = uint8_t;

//! \brief Type for positions within a motif. We expect that motifs are shorter than 65k.
using MotifLen = uint16_t;

//! \brief Type for the score of a motif.
using MotifScore = float;

//! \brief Store {min, mean, max} of a distribution.
struct LengthStat
{
    MotifLen min;
    MotifLen max;
    float mean;
};

//! \brief The boundaries of a stemloop.
typedef std::pair<MotifLen, MotifLen> Coordinate;

//! \brief A loop element in a stemloop.
struct LoopElement
{
    LengthStat length;
    std::vector<profile_char<seqan3::rna4>> profile;
    std::vector<std::unordered_map<MotifLen, SeqNum>> gaps;
    bool is_5prime;
};

//! \brief A stem element in a stemloop.
struct StemElement
{
    LengthStat length;
    std::vector<profile_char<bi_alphabet<seqan3::rna4>>> profile;
    std::vector<std::unordered_map<MotifLen, SeqNum>> gaps;
};

//! \brief A stemloop motif consists of a series of loop and stem elements.
struct StemloopMotif
{
    //! \brief A unique identifier for the motif.
    MotifNum uid;

    //! \brief The position interval, where the stemloop is located in the alignment.
    Coordinate bounds;

    //! \brief The length statistics of the stemloop.
    LengthStat length;

    //! \brief The number of underlying sequences of which the motif was created.
    SeqNum depth;

    //! \brief A vector of loop and stem elements that the stemloop consists of.
    std::vector<std::variant<LoopElement, StemElement>> elements;

    /*!
     * \brief Add a new stem to this motif.
     * \return a reference to the new stem element.
     */
    StemElement & new_stem();

    /*!
     * \brief Add a new loop to this motif.
     * \param is_5prime True for 5' loops, false for 3' loops.
     * \return a reference to the new loop element.
     */
    LoopElement & new_loop(bool is_5prime);

    /*!
     * \brief Analyze the motif's properties based on the MSA and interactions.
     * \param msa The multiple structural alignment.
     */
    void analyze(Msa const & msa);

    /*!
     * \brief Print the motif as RSSP for the Structator program.
     * \param[in,out] os The output stream.
     */
    void print_rssp(std::ofstream & os) const;

    /*!
     * \brief Constructor for a stemloop motif.
     * \param id A unique ID for the motif.
     * \param pos The location of the motif.
     */
    StemloopMotif(MotifNum id, Coordinate pos) :
        uid{id},
        bounds{std::move(pos)},
        length{},
        depth{},
        elements{}
    {}
};

/*!
 * \brief Stream a representation of a stem loop motif.
 * \param os The output stream.
 * \param motif The motif to be printed.
 * \return the stream with the motif representation appended.
 */
std::ostream & operator<<(std::ostream & os, StemloopMotif const & motif);

/*!
 * \brief Create the motif descriptors by analysing a multiple sequence-structure alignment.
 * \param alignment_file The filepath containing the MSA.
 * \param threads The maximum number of threads allowed for execution.
 * \return A vector of motifs.
 */
std::vector<StemloopMotif> create_motifs(std::filesystem::path const & alignment_file, unsigned int threads);

/*!
 * \brief Extract the positions of the stem loops.
 * \param bpseq The base pairing at each position.
 * \param plevel The pseudoknot level at each position.
 * \return a vector of motifs with initialized stemloop positions.
 */
std::vector<StemloopMotif> detect_stemloops(std::vector<int> const & bpseq, std::vector<int> const & plevel);

/*!
 * \brief Analyse whether the profile char is ambiguous.
 * \tparam alph_type The underlying alphabet type of the profile char.
 * \param[in] prof The profile char to check.
 * \return The rank if the profile contains only one entry, the alphabet size otherwise.
 */
template <seqan3::semialphabet alph_type>
unsigned short get_profile_rank(profile_char<alph_type> const & prof)
{
    unsigned short rank = seqan3::alphabet_size<alph_type>;
    for (unsigned short idx = 0; idx < seqan3::alphabet_size<alph_type>; ++idx)
    {
        if (prof.quantity(idx) > 0)
        {
            if (rank == seqan3::alphabet_size<alph_type>)
                rank = idx;
            else
                return seqan3::alphabet_size<alph_type>; // N
        }
    }
    return rank;
}

} // namespace mars
