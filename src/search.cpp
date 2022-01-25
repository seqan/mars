#include <algorithm>
#include <iostream>

#include "search.hpp"
#include "settings.hpp"

namespace mars
{

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
    seqan3::rna4 c = get<0>(stem_item.second);
    bool succ = new_cur.extend_left(c);
    if (succ)
    {
        c = get<1>(stem_item.second);
        succ = new_cur.extend_right(c);
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
    if (scores.size() < settings.xdrop)
        return false;
    else
        return scores.back() < scores[scores.size() - settings.xdrop];
}

void SearchInfo::compute_hits() const
{
    for (auto && [seq, pos] : cursors.back().locate())
        hits.push(seq, pos - motif.bounds.first, motif.uid, scores.back());
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

std::set<MotifLocation, MotifLocationCompare> find_motifs(mars::BiDirectionalIndex const & index,
                                                          std::vector<StemloopMotif> const & motifs)
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

    std::set<MotifLocation, MotifLocationCompare> locations{};
    std::mutex mutex_locations;

    futures.clear();
    for (size_t sidx = 0u; sidx < index.seq_count(); ++sidx)
    {
        if (hits.get(sidx).empty())
            continue;

        futures.push_back(pool->submit([sidx, num_motifs, &mutex_locations, &locations, &hits]
        {
            std::vector<Hit> & hitvec = hits.get(sidx);
            std::sort(hitvec.begin(), hitvec.end(), [] (Hit const & hit1, Hit const & hit2)
            {
                if (hit1.pos != hit2.pos)
                    return hit1.pos < hit2.pos;
                if (hit1.midx != hit2.midx)
                    return hit1.midx < hit2.midx;
                return hit1.score < hit2.score;
            });

            auto left_end = hitvec.cbegin();
            auto right_end = left_end;

            do
            {
                std::vector<float> max_score(num_motifs, 0.f);
                std::vector<Hit> selection{};
                while (right_end != hitvec.cend() && right_end->pos <= left_end->pos + 30)
                {
                    selection.push_back(*right_end);
                    ++right_end;
                }
                std::sort(selection.begin(), selection.end(), [] (Hit const & hit1, Hit const & hit2)
                {
                    if (hit1.midx != hit2.midx)
                        return hit1.midx < hit2.midx;
                    return hit1.score > hit2.score;
                });
                long long base_pos = static_cast<long long>(selection.begin()->pos);
                for (Hit & hit : selection)
                {
                    int pos_diff = static_cast<int>(hit.pos - base_pos);
                    hit.score = static_cast<float>(std::max(0.0, hit.score - (0.1 * pos_diff * pos_diff)));
                    max_score[hit.midx] = std::max(max_score[hit.midx], hit.score);
                }

                uint8_t diversity = 0; // number of different motifs found
                float hit_score = 0;
                for (float score : max_score)
                {
                    diversity += score > 0.f ? 1 : 0;
                    hit_score += score;
                }

                if (diversity > num_motifs / 4 && hit_score * 2 > num_motifs * settings.min_score_per_motif)
                {
                    std::lock_guard<std::mutex> guard(mutex_locations);
                    locations.emplace(hit_score, diversity, base_pos, sidx);
                }

                left_end = right_end;
            } while (right_end != hitvec.cend());
        })); // end of thread pool execution
    }

    for (auto & future : futures)
        future.wait();

    return locations;
}

void print_locations(std::set<MotifLocation, MotifLocationCompare> const & locations,
                     mars::BiDirectionalIndex const & index)
{
    auto print_results = [&locations, &index] (std::ostream & out)
    {
        if (!locations.empty())
            out << " " << std::left << std::setw(35) << "sequence name" << "\t" << "index" << "\t"
                << "pos" << "\t" << "n" << "\t" << "score" << std::endl;
        for (mars::MotifLocation const & loc : locations)
            out << ">" << std::left << std::setw(35) << index.seq_name(loc.sequence) << "\t" << loc.sequence << "\t"
                << loc.position << "\t" << +loc.num_stemloops << "\t" << loc.score << std::endl;
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
