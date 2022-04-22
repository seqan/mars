#pragma once

#include <ostream>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/utility/views/zip.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/types/vector.hpp>
#endif

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

//! \brief Type for sequence id. We expect to have less than 4G sequences.
using SeqNum = size_t;

//! \brief Store {min, mean, max} of a distribution.
struct LengthStat
{
    MotifLen min;
    MotifLen max;
    float mean;

#if SEQAN3_WITH_CEREAL
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(min);
        archive(max);
        archive(mean);
    }
#endif
};

//! \brief The boundaries of a stemloop.
typedef std::pair<MotifLen, MotifLen> Coordinate;

//! \brief A loop element in a stemloop.
struct LoopElement
{
    std::vector<std::vector<std::pair<MotifScore, seqan3::rna4>>> prio;
    std::vector<std::unordered_map<MotifLen, SeqNum>> gaps;
    bool is_5prime;

#if SEQAN3_WITH_CEREAL
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(prio);
        archive(gaps);
        archive(is_5prime);
    }
#endif
};

//! \brief A stem element in a stemloop.
struct StemElement
{
    std::vector<std::vector<std::pair<MotifScore, bi_alphabet<seqan3::rna4>>>> prio;
    std::vector<std::unordered_map<MotifLen, SeqNum>> gaps;

#if SEQAN3_WITH_CEREAL
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(prio);
        archive(gaps);
    }
#endif
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

#if SEQAN3_WITH_CEREAL
    StemloopMotif() : uid{}, bounds{}, length{}, depth{}, elements{}
    {}

    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(uid);
        archive(bounds);
        archive(length);
        archive(depth);
        archive(elements);
    }
#endif
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
 * \param threads The maximum number of threads allowed for execution.
 * \return A vector of motifs.
 */
std::vector<StemloopMotif> create_motifs();

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

/*!
 * \brief Retrieve the log quantities relative to the background distribution.
 * \param pch The profile char.
 * \param depth The number of sequences used to create the profile (for normalization).
 * \return A priority queue with logarithmic scores and the respective RNA characters.
 */
template <seqan3::semialphabet alph_type>
std::vector<std::pair<MotifScore, alph_type>> priority(profile_char<alph_type> pch, size_t depth)
{
    std::vector<std::pair<MotifScore, alph_type>> result{};
    for (auto && [qnt, chr, bg] : seqan3::views::zip(pch.quantities(), pch.alphabet, pch.background_distribution))
    {
        if (qnt > 0)
        {
            result.emplace_back(log2f((qnt + 1.) / pch.one / depth) - bg, chr);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

/*!
 * \brief Write the motifs in rssp format for Structator.
 * \param motifs the motif vector.
 */
void store_rssp(std::vector<StemloopMotif> const & motifs);

#if SEQAN3_WITH_CEREAL
/*!
 * \brief Read motifs from a file.
 * \param motif_file The filename from which where the motifs can be restored.
 * \return the motif vector.
 */
std::vector<StemloopMotif> restore_motifs(std::filesystem::path const & motif_file);

/*!
 * \brief Write the motifs to a file.
 * \param motifs The motifs.
 */
void store_motifs(std::vector<StemloopMotif> const & motifs);
#endif

} // namespace mars
