#include <iomanip>

#include "location.hpp"

namespace mars
{

bool operator<(MotifLocation const & loc1, MotifLocation const & loc2)
{
    if (settings.evalue_filter && loc1.evalue != loc2.evalue)
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

void LocationCollector::push(MotifLocation && loc)
{
    std::lock_guard<std::mutex> guard(mutex_locations);
    emplace_back(loc);
}

void LocationCollector::print()
{
    std::sort(begin(), end());
    if (!settings.result_file.empty())
    {
        logger(1, "Writing " << size() << " results ==> " << settings.result_file << std::endl);
        std::ofstream file_stream(mars::settings.result_file);
        print_results(file_stream);
        file_stream.close();
    }
    else
    {
        logger(1, "Writing " << size() << " results ==> stdout" << std::endl);
        std::lock_guard<std::mutex> guard(mutex_console);
        print_results(std::cout);
    }
}

void LocationCollector::print_results(std::ostream & out)
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
        out << std::left << std::setw(35) << index.seq_name(iter->sequence)
            << "\t" << iter->sequence
            << "\t" << iter->position_start
            << "\t" << iter->position_end
            << "\t" << iter->query_length
            << "\t" << +iter->num_stemloops
            << "\t" << iter->score
            << "\t" << iter->evalue
            << std::endl;
    }
    while (++iter != cend() && (!settings.evalue_filter || iter->evalue < thr));
}

bool operator<(Hit const & hit1, Hit const & hit2)
{
    return hit1.pos < hit2.pos;
}

HitStore::HitStore(size_t seq_count)
{
    hits.resize(seq_count);
}

void HitStore::push(Hit && hit, size_t seq)
{
    std::lock_guard<std::mutex> guard(mutexes[seq % 32]); // distribute the mutexes
    hits[seq].emplace_back(hit);
}

std::vector<Hit> & HitStore::get(size_t seq)
{
    return hits[seq];
}

} // namespace mars
