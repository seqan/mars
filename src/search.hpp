// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <future>
#include <set>
#include <tuple>
#include <variant>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>

#include "index.hpp"
#include "location.hpp"
#include "motif.hpp"

namespace mars
{

//! \brief The type of a stemloop's element iterator.
typedef typename std::vector<std::variant<LoopElement, StemElement>>::const_iterator ElementIter;

//! \brief A storage for futures with concurrent access.
struct ConcurrentFutureVector
{
    std::vector<std::future<void>> futures; //!< The future storage.
    std::mutex mutex; //!< A mutex for concurrent access to `futures`.
};

//! \brief Provides a bi-directional step-by-step stemloop search with backtracking.
class SearchInfo
{
private:
    //! \brief The history of scores and cursors (needed for backtracking)
    std::vector<std::pair<float, seqan3::bi_fm_index_cursor<Index>>> history;

    //! \brief The stemloop to be searched.
    Stemloop const & stemloop;

    //! \brief The resulting stemloop hits are stored here concurrently.
    StemloopHitStore & hits;

    //! \brief Storage for the task futures of locating the hits.
    ConcurrentFutureVector & queries;

public:
    /*!
     * \brief Constructor for a bi-directional search.
     * \param index The index where the search takes place.
     * \param stemloop The stemloop to be searched.
     * \param hits Storage for the resulting stemloop hits.
     * \param queries Storage for the task futures of locating the hits.
     */
    SearchInfo(Index const & index,
               Stemloop const & stemloop,
               StemloopHitStore & hits,
               ConcurrentFutureVector & queries):
        stemloop{stemloop},
        hits{hits},
        queries{queries}
    {
        history.emplace_back(0, index);
    }

    /*!
     * \brief Append a character to the 5' (left) side of the query.
     * \param item The character to be added.
     * \param left Whether the loop is at the 5' side.
     * \returns whether the operation was successful.
     */
    bool append_loop(std::pair<float, seqan3::rna4> item, bool left);

    /*!
     * \brief Append a character pair at both sides of the query.
     * \param stem_item The score and characters to be added.
     * \returns whether the operation was successful.
     */
    bool append_stem(ScoredRnaPair stem_item);

    //! \brief Revert the previous append step, which shrinks the query by one or two characters.
    void backtrack();

    /*!
     * \brief Whether the search should be aborted through the xdrop condition.
     * \return True if the score dropped over the previous x elements, false otherwise.
     */
    [[nodiscard]] bool xdrop() const;

    /*!
     * \brief Determine whether we have reached the last element of the stemloop.
     * \return the end iterator for the stemloop's elements.
     */
    ElementIter stemloop_end() const
    {
        return stemloop.elements.cend();
    }

    //! \brief Locate the current query in the genome and store the result in `hits`.
    void compute_hits() const;
};

/*!
 * \brief Recursive function that descends in the search tree.
 * \tparam MotifElement Type that determines whether we search a loop or stem.
 * \param info A reference to the search information.
 * \param elem_it The current stemloop element where to start the search.
 * \param idx The position in the stemloop element where to start the search.
 */
template <typename MotifElement>
void recurse_search(SearchInfo & info, ElementIter elem_it, Position idx);

/*!
 * \brief Combine hits into motif locations separately for each sequence in range.
 * \param locations The resulting locations.
 * \param hits The hits for each sequence.
 * \param motif The motif, i.e. the vector of stemloops that was subject to the search.
 * \param db_len The total length of all sequences.
 * \param sidx_begin The first sequence in range.
 * \param sidx_end One after the last sequence in range.
 */
void merge_hits(MotifLocationStore & locations,
                StemloopHitStore & hits,
                Motif const & motif,
                size_t db_len,
                size_t sidx_begin,
                size_t sidx_end);

/*!
 * \brief Initiate the recursive search.
 * \param index The index to be searched in.
 * \param motif The motif to be searched.
 */
void find_motif(BiDirectionalIndex const & index, Motif const & motif);

} // namespace mars
