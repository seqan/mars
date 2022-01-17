#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <set>
#include <tuple>
#include <variant>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

#include "motif.hpp"
#include "index.hpp"

namespace mars
{

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

class SearchGenerator
{
private:
    using ElementIter = typename std::vector<std::variant<LoopElement, StemElement>>::const_reverse_iterator;

    BiDirectionalIndex & bds;
    std::vector<std::vector<Hit>> hits;
    std::set<MotifLocation, MotifLocationCompare> locations;

    template <typename MotifElement>
    void recurse_search(StemloopMotif const & motif, ElementIter const & elem_it, MotifLen idx);

public:
    SearchGenerator(BiDirectionalIndex & bds) :
        bds{bds},
        hits{},
        locations{}
    {}

    void find_motifs(std::vector<StemloopMotif> const & motifs, unsigned threads, float min_score);

    std::set<MotifLocation, MotifLocationCompare> const & get_locations() const
    {
        return locations;
    }
};

} // namespace mars
