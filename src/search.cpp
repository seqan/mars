#include <algorithm>
#include <thread>

#ifdef MARS_WITH_OPENMP
    #include <omp.h>
#endif

#include "search.hpp"
#include "settings.hpp"

namespace mars
{

template <typename MotifElement>
void SearchGenerator::recurse_search(StemloopMotif const & motif, ElementIter const & elem_it, MotifLen idx)
{
    if (bds.xdrop())
        return;

    auto const & elem = std::get<MotifElement>(*elem_it);

    if (idx == elem.profile.size())
    {
        auto const next = elem_it + 1;
        if (next == motif.elements.crend())
            bds.compute_hits(hits, motif);
        else if (std::holds_alternative<StemElement>(*next))
            recurse_search<StemElement>(motif, next, 0);
        else
            recurse_search<LoopElement>(motif, next, 0);
        return;
    }

    auto const & prio = elem.prio[elem.profile.size() - idx - 1];

    // try to extend the pattern
    for (auto opt = prio.crbegin(); opt != prio.crend(); ++opt)
    {
        bool succ;
        if constexpr (std::is_same_v<MotifElement, LoopElement>)
            succ = bds.append_loop(*opt, elem.is_5prime);
        else
            succ = bds.append_stem(*opt);

        if (succ)
        {
            recurse_search<MotifElement>(motif, elem_it, idx + 1);
            bds.backtrack();
        }
    }

    // try gaps
    for (auto && [len, num] : elem.gaps[elem.gaps.size() - idx - 1])
        recurse_search<MotifElement>(motif, elem_it, idx + len);
}

void SearchGenerator::find_motifs(std::vector<StemloopMotif> const & motifs, unsigned threads, float min_score)
{
    if (mars::verbose > 0)
        std::cerr << "Start the motif search...\n";
    assert(motifs.size() <= UINT8_MAX);
    uint8_t const num_motifs = motifs.size();
    hits.resize(bds.number_of_seq());

    for (auto const & motif : motifs)
        bds.update_max_offset(motif.bounds.first);

    std::mutex mutex_cerr;

    std::vector<std::thread> thread_pool(1);
    for (StemloopMotif const & motif : motifs)
    {
        thread_pool[0] = std::thread([this, motif, &mutex_cerr]
        {
            if (verbose > 0)
            {
                std::lock_guard<std::mutex> guard(mutex_cerr);
                std::cerr << " Start " << (+motif.uid + 1) << std::endl;
            }
            // start with the hairpin
            auto const iter = motif.elements.crbegin();
            if (std::holds_alternative<LoopElement>(*iter))
                recurse_search<LoopElement>(motif, iter, 0);
            else
                recurse_search<StemElement>(motif, iter, 0);
            if (verbose > 0)
            {
                std::lock_guard<std::mutex> guard(mutex_cerr);
                std::cerr << " Finish " << (+motif.uid + 1) << std::endl;
            }
        });
        thread_pool[0].join();
    }

//    for (auto & thr : thread_pool)
//        thr.join();

    std::for_each(hits.begin(), hits.end(), [] (std::vector<Hit> & hitvec)
    {
        std::sort(hitvec.begin(), hitvec.end(), [] (Hit const & a, Hit const & b)
        {
            if (a.pos != b.pos)
                return a.pos < b.pos;
            if (a.midx != b.midx)
                return a.midx < b.midx;
            return a.score < b.score;
        });
    });

    size_t num_results = 0;
    #pragma omp parallel for num_threads(threads)
    for (size_t sidx = 0u; sidx < bds.number_of_seq(); ++sidx)
    {
        std::vector<Hit> const & hitvec = hits[sidx];
        if (hitvec.empty())
            continue;

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
            std::sort(selection.begin(), selection.end(), [] (Hit const & a, Hit const & b)
            {
                if (a.midx != b.midx)
                    return a.midx < b.midx;
                return a.score > b.score;
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

            base_pos -= static_cast<long long>(bds.get_max_offset());

            if (diversity > 0 && hit_score > num_motifs * min_score)
            {
                #pragma omp critical (nres)
                {
                    locations.emplace(hit_score, diversity, base_pos, sidx);
                    ++num_results;
                }
            }

            left_end = right_end;
        } while (right_end != hitvec.cend());
    }
    if (verbose > 1)
        std::cerr << "Found " << num_results << " matches." << std::endl;
}

} // namespace mars
