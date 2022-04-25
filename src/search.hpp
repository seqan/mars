#pragma once

#include <array>
#include <cmath>
#include <cstdint>
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

using ElementIter = typename std::vector<std::variant<LoopElement, StemElement>>::const_iterator;

void find_motifs(BiDirectionalIndex const & index, std::vector<StemloopMotif> const & motifs);

void merge_hits(LocationCollector & locations,
                HitStore & hits,
                std::vector<StemloopMotif> const & motifs,
                size_t db_len,
                size_t sidx_begin,
                size_t sidx_end);

typedef std::pair<float, seqan3::bi_fm_index_cursor<Index>> ScoredCursor;

//! \brief Provides a bi-directional search step-by-step with backtracking.
struct SearchInfo
{
private:
    StemloopMotif const & motif;
    HitStore & hits;
    std::vector<std::future<void>> & futures;

    //! \brief The history of scores and cursors (needed for backtracking)
    std::vector<ScoredCursor> history;

public:
    /*!
     * \brief Constructor for a bi-directional search.
     */
    SearchInfo(Index const & index, StemloopMotif const & stemloop, HitStore & hits,
               std::vector<std::future<void>> & futures):
        motif{stemloop},
        hits{hits},
        futures{futures},
        history{}
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
    bool append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item);

    //! \brief Revert the previous append step, which shrinks the query by one or two characters.
    void backtrack();

    /*!
     * \brief Whether the search should be aborted through the xdrop condition.
     * \return True if the score dropped over the previous x elements, false otherwise.
     */
    [[nodiscard]] bool xdrop() const;

    ElementIter motif_end() const
    {
        return motif.elements.cend();
    }

    /*!
     * \brief Perform the search with the current query and store the result in `hits`.
     */
    void compute_hits() const;
};

template <typename MotifElement>
void recurse_search(SearchInfo & info, ElementIter const elem_it, MotifLen idx);

} // namespace mars
