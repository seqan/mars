// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <ostream>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

#include <seqan3/alphabet/gap/gapped.hpp>
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

//! \brief Type for positions within a stemloop. We expect that stemloops are shorter than 65k.
using Position = uint16_t;

//! \brief The boundaries of a stemloop.
typedef std::pair<Position, Position> Bounds;

//! \brief A pair of score and gapped RNA bi-character (to represent stems).
typedef std::pair<float, bi_alphabet<seqan3::gapped<seqan3::rna4>>> ScoredRnaPair;

//! \brief A pair of score and RNA character (to represent loops).
typedef std::pair<float, seqan3::rna4> ScoredRna;

//! \brief A loop element in a stemloop.
struct LoopElement
{
    //! \brief Prioritized list of loop characters.
    std::vector<std::vector<ScoredRna>> prio;

    //! \brief For each position we store the length and number of gaps in a hash map.
    std::vector<std::unordered_map<Position, size_t>> gaps;

    //! \brief Flag whether the loop is on the left side (towards 5').
    bool leftsided;

#if SEQAN3_WITH_CEREAL
    /*!
     * \brief Function that guides Cereal serialization.
     * \tparam Archive Type of the Cereal archive.
     * \param archive The archive.
     */
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(prio);
        archive(gaps);
        archive(leftsided);
    }
#endif
};

//! \brief A stem element in a stemloop.
struct StemElement
{
    //! \brief Prioritized list of stem characters.
    std::vector<std::vector<ScoredRnaPair>> prio;

    //! \brief For each position we store the length and number of gaps in a hash map.
    std::vector<std::unordered_map<Position, size_t>> gaps;

#if SEQAN3_WITH_CEREAL
    /*!
     * \brief Function that guides Cereal serialization.
     * \tparam Archive Type of the Cereal archive.
     * \param archive The archive.
     */
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(prio);
        archive(gaps);
    }
#endif
};

//! \brief A stemloop consists of a series of loop and stem elements.
struct Stemloop
{
    //! \brief A unique identifier for the stemloop.
    uint8_t uid;

    //! \brief The position interval, where the stemloop is located in the alignment.
    Bounds bounds;

    //! \brief The minimum and maximum length of the stemloop.
    Bounds length;

    //! \brief A vector of loop and stem elements that the stemloop consists of.
    std::vector<std::variant<LoopElement, StemElement>> elements;

    /*!
     * \brief Constructor for a stemloop.
     * \param id A unique ID for the stemloop.
     * \param pos The location of the stemloop.
     */
    Stemloop(uint8_t id, Bounds pos) :
        uid{id},
        bounds{std::move(pos)},
        length{},
        elements{}
    {}

#if SEQAN3_WITH_CEREAL
    //! \brief Default constructor for serialization.
    Stemloop() : uid{}, bounds{}, length{}, elements{}
    {}

    /*!
     * \brief Function that guides Cereal serialization.
     * \tparam Archive Type of the Cereal archive.
     * \param archive The archive.
     */
    template <seqan3::cereal_archive Archive>
    void serialize(Archive & archive)
    {
        archive(uid);
        archive(bounds);
        archive(length);
        archive(elements);
    }
#endif

    /*!
     * \brief Analyze the stemloop's properties based on the MSA and interactions.
     * \param msa The multiple structural alignment.
     */
    void analyze(Msa const & msa);

    /*!
     * \brief Print the stemloop as RSSP for the Structator program.
     * \param[in,out] os The output stream.
     */
    void print_rssp(std::ofstream & os) const;
};

//! \brief A motif is a collection of stemloops.
typedef std::vector<Stemloop> Motif;

/*!
 * \brief Stream a representation of a stem loop.
 * \param os The output stream.
 * \param stemloop The stemloop to be printed.
 * \return the stream with the stemloop representation appended.
 */
std::ostream & operator<<(std::ostream & os, Stemloop const & stemloop);

/*!
 * \brief Write the motif in rssp format for Structator.
 * \param motif The motif.
 */
void store_rssp(Motif const & motif);

/*!
 * \brief Create the motif descriptors by analysing a multiple sequence-structure alignment.
 * \param threads The maximum number of threads allowed for execution.
 * \return A motif (vector of stemloops).
 */
Motif create_motif();

/*!
 * \brief Extract the positions of the stem loops.
 * \param bpseq The base pairing at each position.
 * \param plevel The pseudoknot level at each position.
 * \return a motif with initialized stemloop positions.
 */
Motif detect_stemloops(std::vector<int> const & bpseq, std::vector<int> const & plevel);

/*!
 * \brief Retrieve the log quantities relative to the background distribution.
 * \param pch The profile char.
 * \param depth The number of sequences used to create the profile (for normalization).
 * \return A priority queue with logarithmic scores and the respective RNA characters.
 */
template <seqan3::semialphabet alph_type>
std::vector<std::pair<float, alph_type>> priority(profile_char<alph_type> pch, size_t depth)
{
    std::vector<std::pair<float, alph_type>> result{};
    for (auto && [qnt, chr, bg] : seqan3::views::zip(pch.quantities(), pch.alphabet, pch.background_distribution))
        if (qnt > 0)
            result.emplace_back(log2f((qnt + 1.) / pch.one / depth) - bg, chr);
    std::sort(result.begin(), result.end());
    return result;
}

#if SEQAN3_WITH_CEREAL
/*!
 * \brief Read motif from a file.
 * \param motif_file The filename from which where the motif can be restored.
 * \return the motif.
 */
Motif restore_motif(std::filesystem::path const & motif_file);

/*!
 * \brief Write the motif to a file.
 * \param motifs The motif.
 */
void store_motif(Motif const & motif);
#endif

} // namespace mars
