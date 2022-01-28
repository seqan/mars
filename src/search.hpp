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
    double evalue;
    float score;
    uint8_t num_stemloops;
    size_t position_start;
    size_t position_end;
    size_t query_length;
    size_t sequence;

    MotifLocation(double eval, float sco, uint8_t num, size_t spos, size_t epos, size_t len, size_t seq):
        evalue{eval}, score{sco}, num_stemloops{num}, position_start{spos}, position_end{epos}, query_length{len},
        sequence{seq}
    {}
};

struct MotifLocationCompare {
    bool operator()(MotifLocation const & loc1, MotifLocation const & loc2) const
    {
        if (loc1.evalue != loc2.evalue)
            return loc1.evalue < loc2.evalue;
        if (loc1.score != loc2.score)
            return loc1.score > loc2.score;
        if (loc1.num_stemloops != loc2.num_stemloops)
            return loc1.num_stemloops > loc2.num_stemloops;
        if (loc1.query_length != loc2.query_length)
            return loc1.query_length > loc2.query_length;
        if (loc1.sequence != loc2.sequence)
            return loc1.sequence < loc2.sequence;
        if (loc1.position_start != loc2.position_start)
            return loc1.position_start < loc2.position_start;
        return loc1.position_end < loc2.position_end;
    }
};

using SortedLocations = std::set<MotifLocation, MotifLocationCompare>;

struct Hit
{
    long long pos;
    long long length;
    uint8_t midx;
    float score;

    Hit(long long pos, long long length, uint8_t midx, float score) : pos{pos}, length{length}, midx{midx}, score{score} {}
};

bool operator<(Hit const & hit1, Hit const & hit2);

std::ostream & operator<<(std::ostream & ostr, Hit const & hit);

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

    void push(size_t seq, long long pos, long long length, uint8_t midx, float score)
    {
        std::lock_guard<std::mutex> guard(mutex_hits);
        hits[seq].emplace_back(pos, length, midx, score);
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
