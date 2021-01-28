#pragma once

#include "motif.hpp"
#include "index.hpp"

namespace mars
{

struct hit
{
    uint8_t id;
    size_t seq;
    size_t pos;
    float score;
};

class search_generator
{
private:
    bi_directional_search bds;
    std::vector<hit> hits;

public:
    explicit search_generator(index_type index) : bds{std::move(index)}, hits{}
    {}

    void find_motif(stemloop_motif const & motif);
};

} // namespace mars
