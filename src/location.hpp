// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <ostream>
#include <vector>

#include "index.hpp"
#include "settings.hpp"

namespace mars
{

struct MotifLocation
{
    double evalue;
    float score;
    uint8_t num_stemloops;
    size_t position_start;
    size_t position_end;
    size_t query_length;
    size_t sequence;
};

bool operator<(MotifLocation const & loc1, MotifLocation const & loc2);

class LocationCollector : public std::vector<MotifLocation>
{
public:
    explicit LocationCollector(BiDirectionalIndex const & index) : index{index} {}
    void print();
    void push(MotifLocation && loc);
private:
    BiDirectionalIndex const & index;
    void print_results(std::ostream & out);
    std::mutex mutex_locations;
};

struct Hit
{
    long long pos;
    long long length;
    uint8_t midx;
    float score;
};

bool operator<(Hit const & hit1, Hit const & hit2);

struct HitStore
{
private:
    std::vector<std::vector<Hit>> hits;
    std::array<std::mutex, 256> mutexes;

public:
    explicit HitStore(size_t seq_count);
    void push(Hit && hit, size_t seq);
    std::vector<Hit> & get(size_t seq);
};

} // namespace mars
