// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <iostream>

#include <seqan3/utility/parallel/detail/latch.hpp>

#include "search.hpp"
#include "settings.hpp"

namespace mars
{

bool SearchInfo::append_loop(std::pair<float, seqan3::rna4> item, bool left)
{
    bool succ;
    seqan3::bi_fm_index_cursor<Index> new_cur(history.back().second);

    if (left)
        succ = new_cur.extend_left(item.second);
    else
        succ = new_cur.extend_right(item.second);

    if (succ)
        history.emplace_back(history.back().first + item.first, new_cur);
    return succ;
}

bool SearchInfo::append_stem(ScoredRnaPair stem_item)
{
    seqan3::bi_fm_index_cursor<Index> new_cur(history.back().second);
    bool succ = true;
    if (stem_item.second.first() != seqan3::gap())
        succ = new_cur.extend_left(stem_item.second.first().convert_unsafely_to<seqan3::rna4>());
    if (succ && stem_item.second.second() != seqan3::gap())
        succ = new_cur.extend_right(stem_item.second.second().convert_unsafely_to<seqan3::rna4>());
    if (succ)
        history.emplace_back(history.back().first + stem_item.first, new_cur);
    return succ;
}

void SearchInfo::backtrack()
{
    history.pop_back();
}

bool SearchInfo::xdrop() const
{
    if (history.back().second.query_length() > stemloop.length.second)
        return true;
    if (history.size() < settings.xdrop)
        return false;
    else
        return history.back().first < history[history.size() - settings.xdrop].first;
}

void SearchInfo::compute_hits() const
{
    auto score = history.back().first;
    auto cur = history.back().second;
    auto const len = static_cast<long long>(cur.query_length());
    if (len >= stemloop.length.first && len > 5 && score > 0)
    {
        std::lock_guard<std::mutex> guard(queries.mutex);
        queries.futures.push_back(pool->submit([cur, store = &hits, off = stemloop.bounds.first, len,
                                                uid = stemloop.uid, score]
        {
            for (auto && [seq, pos] : cur.locate())
                store->push({static_cast<long long>(pos) - off, len, uid, score}, seq);
        }));
    }
}

template <typename MotifElement>
void recurse_search(SearchInfo & info, ElementIter elem_it, Position idx)
{
    if (info.xdrop())
        return;

    auto const & elem = std::get<MotifElement>(*elem_it);

    if (idx == elem.prio.size())
    {
        auto const next = elem_it + 1;
        if (next == info.stemloop_end())
            info.compute_hits();
        else if (std::holds_alternative<StemElement>(*next))
            recurse_search<StemElement>(info, next, 0);
        else
            recurse_search<LoopElement>(info, next, 0);
        return;
    }

    auto const & prio = elem.prio[idx];

    // try to extend the pattern
    for (auto opt = prio.crbegin(); opt != prio.crend(); ++opt)
    {
        bool succ;
        if constexpr (std::is_same_v<MotifElement, LoopElement>)
            succ = info.append_loop(*opt, elem.leftsided);
        else
            succ = info.append_stem(*opt);

        if (succ)
        {
            recurse_search<MotifElement>(info, elem_it, idx + 1);
            info.backtrack();
        }
    }

    // try gaps
    for (auto const & len_num : elem.gaps[idx])
        recurse_search<MotifElement>(info, elem_it, idx + len_num.first);
}

void find_motif(mars::BiDirectionalIndex const & index, Motif const & motif)
{
    StemloopHitStore hits(index.get_names().size());

    logger(1, "Stem loop search...");
    assert(motif.size() <= UINT8_MAX);
    uint8_t const num_motifs = motif.size();

    ConcurrentFutureVector queries;
    std::vector<std::future<void>> search_tasks;
    seqan3::detail::latch lat{num_motifs};
    for (size_t idx = 0; idx < num_motifs; ++idx)
    {
        search_tasks.push_back(pool->submit([&index, &motif, &hits, &queries, &lat, idx]
        {
            // initiate recursive search
            SearchInfo info(index.raw(), motif[idx], hits, queries);
            auto const iter = motif[idx].elements.cbegin();
            lat.wait();
            if (std::holds_alternative<LoopElement>(*iter))
                recurse_search<LoopElement>(info, iter, 0);
            else
                recurse_search<StemElement>(info, iter, 0);
            logger(1, " " << (idx + 1));
        }));
        lat.arrive();
    }
    for (auto & future : search_tasks)
        future.wait();
    logger(1, "\nWaiting for " << queries.futures.size() << " queries to complete...");
    std::chrono::steady_clock::time_point tm0 = std::chrono::steady_clock::now();
    for (auto & future : queries.futures)
        future.wait();
    queries.futures.clear();
    auto const sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - tm0).count();
    logger(1, " finished (" << sec << "s)." << std::endl);

    // collect the hits asynchronously and print
    MotifLocationStore locations(index.get_names());
    size_t const seqnum = index.get_names().size();
    size_t const db_len = index.raw().size() - (seqnum > 1 ? seqnum : 2);
    size_t const delta = (seqnum - 1) / settings.nthreads + 1; // ceil
    std::vector<std::future<void>> futures;
    for (size_t sidx = 0; sidx < seqnum; sidx += delta)
    {
        futures.push_back(pool->submit(merge_hits, std::ref(locations), std::ref(hits), std::ref(motif),
                                       db_len, sidx, std::min(sidx + delta, seqnum)));
    }
    for (auto & future : futures)
        future.wait();

    locations.print();
}

void merge_hits(MotifLocationStore & locations,
                StemloopHitStore & hits,
                Motif const & motif,
                size_t db_len,
                size_t sidx_begin,
                size_t sidx_end)
{
    for (size_t sidx = sidx_begin; sidx < sidx_end; ++sidx)
    {
        std::vector<StemloopHit> & hitvec = hits.get(sidx);
        if (hitvec.empty())
            continue;

        std::sort(hitvec.begin(), hitvec.end()); // sort by genome position
        auto left_end = hitvec.cbegin();
        auto right_end = left_end;
        auto const stop = hitvec.cend();

        do
        {
            std::vector<std::vector<StemloopHit>::const_iterator> best_hits(motif.size(), stop);
            // we allow a position divergence of half alignment length
            while (right_end != stop && right_end->pos <= left_end->pos + motif.back().bounds.second / 2)
            {
                auto & iter = best_hits[right_end->midx];
                if (iter == stop || iter->score < right_end->score)
                    iter = right_end;
                ++right_end;
            }

            size_t pos_min{LLONG_MAX};
            size_t pos_max{0};
            size_t query_len{0};
            float bit_score{0};
            uint8_t diversity{0}; // number of different stemloops found

            for (auto & hit : best_hits)
            {
                if (hit != stop)
                {
                    pos_min = std::min(static_cast<size_t>(hit->pos + motif[hit->midx].bounds.first), pos_min);
                    pos_max = std::max(static_cast<size_t>(hit->pos + hit->length + motif[hit->midx].bounds.first),
                                       pos_max);
                    bit_score += hit->score;
                    query_len += hit->length;
                    ++diversity;
                }
            }

            if (std::isnan(settings.score_filter) || // evalue filter
                (diversity > motif.size() / 4 &&
                 bit_score > static_cast<float>(motif.size()) * settings.score_filter))
            {
                double const evalue = static_cast<double>(db_len * query_len) / exp2(bit_score);
                locations.push({evalue, bit_score, diversity, pos_min, pos_max, query_len, sidx});
            }

            left_end = right_end;
        } while (right_end != stop);
    }
}

} // namespace mars
