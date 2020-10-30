#pragma once

#include <vector>

namespace mars
{

//!\brief Store the type of statistics values.
typedef uint32_t value_type;

//!\brief Store {start, end} positions of an inverval.
typedef std::pair<value_type, value_type> interval_type;

//!\brief Store {min, median, max} of a distribution.
typedef std::tuple<value_type, value_type, value_type> stat_type;

class secondary_structure
{
public:
    //enum loop_type {stem, hairpin, bulge, internal};
    std::string loop_fwd;
    std::string loop_bwd;
    std::string stem;

    stat_type length;
    // TODO: gaps
};

class stem_loop_motif
{
private:
    interval_type bounds;
    stat_type length;
    std::vector<secondary_structure> elements;
};

class motif_store
{
private:
    std::vector<int> bpseq;
    std::vector<std::pair<int, int>> stemloops;

public:
    explicit motif_store(std::vector<int> bpseq): bpseq(std::move(bpseq))
    {}

    void stem_loop_partition(std::vector<int> plevel);

    std::vector<std::pair<int, int>> get_stemloops()
    {
        return stemloops;
    }
};

} // namespace mars
