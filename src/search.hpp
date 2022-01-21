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

#include "motif.hpp"
#include "index.hpp"

namespace mars
{

using ElementIter = typename std::vector<std::variant<LoopElement, StemElement>>::const_reverse_iterator;

struct MotifLocation
{
    float score;
    uint8_t num_stemloops;
    long long position;
    size_t sequence;

    MotifLocation(float s, uint8_t n, long long p, size_t i):
        score{s}, num_stemloops{n}, position{p}, sequence{i}
    {}
};

struct MotifLocationCompare {
    bool operator()(MotifLocation const & a, MotifLocation const & b) const
    {
        if (a.score != b.score)
            return a.score > b.score;
        if (a.num_stemloops != b.num_stemloops)
            return a.num_stemloops > b.num_stemloops;
        if (a.sequence != b.sequence)
            return a.sequence < b.sequence;
        return a.position < b.position;
    }
};

struct Hit
{
    long long pos;
    uint8_t midx;
    float score;

    Hit(long long pos, uint8_t midx, float score) : pos{pos}, midx{midx}, score{score} {}
};

struct HitStore
{
private:
    std::vector<std::vector<Hit>> hits;
    std::mutex mutex_hits;

public:
    explicit HitStore(size_t seq_count)
    {
        hits.resize(seq_count);
    }

    void push(size_t seq, long long pos, uint8_t midx, float score)
    {
        std::lock_guard<std::mutex> guard(mutex_hits);
        hits[seq].emplace_back(pos, midx, score);
    }

    std::vector<Hit> & get(size_t idx)
    {
        return hits[idx];
    }
};

std::set<MotifLocation, MotifLocationCompare> find_motifs(mars::BiDirectionalIndex const & index,
                                                          std::vector<StemloopMotif> const & motifs);

//! \brief Provides a bi-directional search step-by-step with backtracking.
struct SearchInfo
{
private:
    StemloopMotif const & motif;
    HitStore & hits;

    //! \brief The history of cursors (needed for backtracking).
    std::vector<seqan3::bi_fm_index_cursor<Index>> cursors;

    //! \brief The history of scores;
    std::vector<float> scores;

public:
    /*!
     * \brief Constructor for a bi-directional search.
     */
    SearchInfo(Index const & index, StemloopMotif const & stemloop, HitStore & hits):
        motif{stemloop},
        hits{hits},
        cursors{},
        scores{}
    {
        cursors.emplace_back(index);
        scores.emplace_back(0);
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
        return motif.elements.crend();
    }

    /*!
     * \brief Perform the search with the current query and store the result in `hits`.
     */
    void compute_hits() const;
};

template <typename MotifElement>
void recurse_search(SearchInfo & info, ElementIter const elem_it, MotifLen idx);

void print_locations(std::set<MotifLocation, MotifLocationCompare> const & locations,
                     mars::BiDirectionalIndex const & index);

} // namespace mars
