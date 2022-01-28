#include <algorithm>
#include <iostream>

#include "search.hpp"
#include "settings.hpp"

namespace mars
{

std::ostream & operator<<(std::ostream & ostr, Hit const & hit)
{
    return ostr << "(" << +hit.midx << "|" << hit.pos << ", " << hit.length << ", " << hit.score << ")";
}

bool operator<(Hit const & hit1, Hit const & hit2)
{
    return hit1.pos < hit2.pos;
}

bool SearchInfo::append_loop(std::pair<float, seqan3::rna4> item, bool left)
{
    bool succ;
    seqan3::bi_fm_index_cursor<Index> new_cur(cursors.back());

    if (left)
        succ = new_cur.extend_left(item.second);
    else
        succ = new_cur.extend_right(item.second);

    if (succ)
    {
        cursors.push_back(new_cur);
        scores.push_back(scores.back() + item.first);
    }
    return succ;
}

bool SearchInfo::append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item)
{
    seqan3::bi_fm_index_cursor<Index> new_cur(cursors.back());
    using seqan3::get;
    seqan3::rna4 chr = get<0>(stem_item.second);
    bool succ = new_cur.extend_left(chr);
    if (succ)
    {
        chr = get<1>(stem_item.second);
        succ = new_cur.extend_right(chr);
    }
    if (succ)
    {
        cursors.push_back(new_cur);
        scores.push_back(scores.back() + stem_item.first);
    }
    return succ;
}

void SearchInfo::backtrack()
{
    scores.pop_back();
    cursors.pop_back();
}

bool SearchInfo::xdrop() const
{
    if (cursors.back().query_length() > motif.length.max)
        return true;
    if (scores.size() < settings.xdrop)
        return false;
    else
        return scores.back() < scores[scores.size() - settings.xdrop];
}

void SearchInfo::compute_hits() const
{
    long long const len = static_cast<long long>(cursors.back().query_length());
    if (len > 0 && scores.back() > 0)
        for (auto && [seq, pos] : cursors.back().locate())
            hits.push(seq, static_cast<long long>(pos) - motif.bounds.first, len, motif.uid, scores.back());
}

template <typename MotifElement>
void recurse_search(SearchInfo & info, ElementIter const elem_it, MotifLen idx)
{
    if (info.xdrop())
        return;

    auto const & elem = std::get<MotifElement>(*elem_it);

    if (idx == elem.profile.size())
    {
        auto const next = elem_it + 1;
        if (next == info.motif_end())
            info.compute_hits();
        else if (std::holds_alternative<StemElement>(*next))
            recurse_search<StemElement>(info, next, 0);
        else
            recurse_search<LoopElement>(info, next, 0);
        return;
    }

    auto const & prio = elem.prio[elem.profile.size() - idx - 1];

    // try to extend the pattern
    for (auto opt = prio.crbegin(); opt != prio.crend(); ++opt)
    {
        bool succ;
        if constexpr (std::is_same_v<MotifElement, LoopElement>)
            succ = info.append_loop(*opt, elem.is_5prime);
        else
            succ = info.append_stem(*opt);

        if (succ)
        {
            recurse_search<MotifElement>(info, elem_it, idx + 1);
            info.backtrack();
        }
    }

    // try gaps
    for (auto const & len_num : elem.gaps[elem.gaps.size() - idx - 1])
        recurse_search<MotifElement>(info, elem_it, idx + len_num.first);
}

SortedLocations find_motifs(mars::BiDirectionalIndex const & index, std::vector<StemloopMotif> const & motifs)
{
    HitStore hits(index.seq_count());

    logger(1, "Stem loop search...");
    assert(motifs.size() <= UINT8_MAX);
    uint8_t const num_motifs = motifs.size();

    std::vector<std::future<void>> futures;
    for (size_t idx = 0; idx < num_motifs; ++idx)
    {
        futures.push_back(pool->submit([idx, &index, &motifs, &hits]
        {
            SearchInfo info(index.raw(), motifs[idx], hits);

            // start with the hairpin
            auto const iter = motifs[idx].elements.crbegin();
            if (std::holds_alternative<LoopElement>(*iter))
                recurse_search<LoopElement>(info, iter, 0);
            else
                recurse_search<StemElement>(info, iter, 0);
            logger(1, " " << (idx + 1));
        }));
    }
    for (auto & future : futures)
        future.wait();
    logger(1, " finished." << std::endl);

    SortedLocations locations{};
    std::mutex mutex_locations;
    size_t const db_len = index.raw().size() - (index.seq_count() > 1 ? index.seq_count() : 2);

    futures.clear();
    for (size_t sidx = 0u; sidx < index.seq_count(); ++sidx)
    {
        if (hits.get(sidx).empty())
            continue;

        futures.push_back(pool->submit([sidx, &motifs, &mutex_locations, &locations, &hits, db_len]
        {
            std::vector<Hit> & hitvec = hits.get(sidx);
            std::sort(hitvec.begin(), hitvec.end()); // sort by genome position
            auto left_end = hitvec.cbegin();
            auto right_end = left_end;
            auto const stop = hitvec.cend();

            do
            {
                std::vector<std::vector<Hit>::const_iterator> best_hits(motifs.size(), stop);
                // we allow a position divergence of half alignment length
                while (right_end != stop && right_end->pos <= left_end->pos + motifs.back().bounds.second / 2)
                {
                    auto & iter = best_hits[right_end->midx];
                    if (iter == stop || iter->score < right_end->score)
                        iter = right_end;
                    ++right_end;
                }

                long long pos_min{LLONG_MAX};
                long long pos_max{0};
                long long query_len{0};
                float bit_score{0};
                uint8_t diversity{0}; // number of different motifs found

                for (auto & hit : best_hits)
                {
                    if (hit != stop)
                    {
                        pos_min = std::min(hit->pos + motifs[hit->midx].bounds.first, pos_min);
                        pos_max = std::max(hit->pos + hit->length + motifs[hit->midx].bounds.first, pos_max);
                        bit_score += hit->score;
                        query_len += hit->length;
                        ++diversity;
                    }
                }

                double evalue = static_cast<double>(db_len * query_len) / exp2(bit_score);
                if (evalue <= settings.max_evalue)
                {
                    std::lock_guard<std::mutex> grd(mutex_locations);
                    locations.emplace(evalue, bit_score, diversity, pos_min, pos_max, query_len, sidx);
                }

                left_end = right_end;
            } while (right_end != stop);
        })); // end of thread pool execution
    }

    for (auto & future : futures)
        future.wait();

    return locations;
}

void print_locations(SortedLocations const & locations, BiDirectionalIndex const & index)
{
    auto print_results = [&locations, &index] (std::ostream & out)
    {
        if (!locations.empty())
            out << std::left << std::setw(35) << "sequence name" << "\t" << "index" << "\t"
                << "pos" << "\t" << "end" << "\t" << "qlen" << "\t" << "n" << "\t" << "score" << "\t" << "e-value"
                << std::endl;
        for (mars::MotifLocation const & loc : locations)
            out << std::left << std::setw(35) << index.seq_name(loc.sequence) << "\t" << loc.sequence << "\t"
                << loc.position_start << "\t" << loc.position_end << "\t" << loc.query_length << "\t"
                << +loc.num_stemloops << "\t" << loc.score << "\t" << loc.evalue << std::endl;
    };

    if (!settings.result_file.empty())
    {
        logger(1, "Writing " << locations.size() << " results ==> " << mars::settings.result_file << std::endl);
        std::ofstream file_stream(mars::settings.result_file);
        print_results(file_stream);
        file_stream.close();
    }
    else
    {
        logger(1, "Writing " << locations.size() << " results ==> stdout" << std::endl);
        print_results(std::cout);
    }
}

} // namespace mars
