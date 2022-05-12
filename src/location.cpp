// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#include <iomanip>
#include <fstream>

#include "location.hpp"
#include "settings.hpp"

namespace mars
{

bool operator<(MotifLocation const & lhs, MotifLocation const & rhs)
{
    if (std::isnan(settings.score_filter) && lhs.evalue != rhs.evalue)
        return lhs.evalue < rhs.evalue;
    if (lhs.score != rhs.score)
        return lhs.score > rhs.score;
    if (lhs.num_stemloops != rhs.num_stemloops)
        return lhs.num_stemloops > rhs.num_stemloops;
    if (lhs.query_length != rhs.query_length)
        return lhs.query_length > rhs.query_length;
    if (lhs.sequence != rhs.sequence)
        return lhs.sequence < rhs.sequence;
    if (lhs.position_start != rhs.position_start)
        return lhs.position_start < rhs.position_start;
    return lhs.position_end < rhs.position_end;
}

void MotifLocationStore::push(MotifLocation && loc)
{
    std::lock_guard<std::mutex> guard(mutex_locations);
    emplace_back(loc);
}

void MotifLocationStore::print(std::ostream & out)
{
    out << std::left << std::setw(35) << "sequence name"
        << "\t" << "index"
        << "\t" << "pos"
        << "\t" << "end"
        << "\t" << "qlen"
        << "\t" << "n"
        << "\t" << "score"
        << "\t" << "e-value"
        << std::endl;

    std::lock_guard<std::mutex> guard(mutex_locations);
    if (empty())
        return;

    auto iter = cbegin();
    double const thr = std::max(std::sqrt(iter->evalue) * 10, 1e-10);
    do
    {
        out << std::left << std::setw(35) << names[iter->sequence]
            << "\t" << iter->sequence
            << "\t" << iter->position_start
            << "\t" << iter->position_end
            << "\t" << iter->query_length
            << "\t" << +iter->num_stemloops
            << "\t" << iter->score
            << "\t" << iter->evalue
            << std::endl;
    }
    while (++iter != cend() && (!std::isnan(settings.score_filter) || iter->evalue < thr));
}

void MotifLocationStore::print()
{
    std::sort(begin(), end());
    if (!settings.result_file.empty())
    {
        logger(1, "Writing the best of " << size() << " results ==> " << settings.result_file << std::endl);
        std::ofstream file_stream(mars::settings.result_file);
        print(file_stream);
        file_stream.close();
    }
    else
    {
        logger(1, "Writing the best of " << size() << " results ==> stdout" << std::endl);
        std::lock_guard<std::mutex> guard(mutex_console);
        print(std::cout);
    }
}

bool operator<(StemloopHit const & lhs, StemloopHit const & rhs)
{
    return lhs.pos < rhs.pos;
}

void StemloopHitStore::push(StemloopHit && hit, size_t seq)
{
    std::lock_guard<std::mutex> guard(mutexes[seq % 256]); // distribute the mutexes
    hits[seq].emplace_back(hit);
}

std::vector<StemloopHit> & StemloopHitStore::get(size_t seq)
{
    return hits[seq];
}

} // namespace mars
