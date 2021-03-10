#include <seqan3/range/views/zip.hpp>

#ifdef MARS_WITH_OPENMP
    #include <omp.h>
#endif

#include "search.hpp"

namespace mars
{

template <seqan3::semialphabet Alphabet>
inline std::set<std::pair<MotifScore, Alphabet>> SearchGenerator::priority(profile_char<Alphabet> const & prof) const
{
    std::set<std::pair<MotifScore, Alphabet>> result{};

    auto const & quantities = prof.log_quantities();
    auto const & bg = background_distr.get<seqan3::alphabet_size<Alphabet>>();

    for (auto && [idx, score, bg] : seqan3::views::zip(std::ranges::views::iota(0), quantities, bg))
        if (score - bg - log_depth >= -2.f)
            result.emplace(score - bg - log_depth, Alphabet{}.assign_rank(idx));

    return std::move(result);
}

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

    auto prio = priority(elem.profile[elem.profile.size() - idx - 1]);

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

void SearchGenerator::find_motifs(std::vector<StemloopMotif> const & motifs)
{
    assert(motifs.size() <= UINT8_MAX);
    uint8_t const num_motifs = motifs.size();
    hits.resize(bds.get_num_seq());

//    #pragma omp parallel for num_threads(2)
    for (uint8_t midx = 0; midx < num_motifs; ++midx)
    {
        auto const & motif = motifs[midx];
        // start with the hairpin
        auto const iter = motif.elements.crbegin();
        recurse_search<LoopElement>(motif, iter, 0);
        std::cerr << motif << std::endl;
    }

    #pragma omp parallel for num_threads(4)
    for (uint16_t sidx = 0u; sidx < bds.get_num_seq(); ++sidx)
        std::sort(hits[sidx].begin(), hits[sidx].end(), [] (Hit const & a, Hit const & b)
        {
            if (a.pos != b.pos)
                return a.pos < b.pos;
            if (a.midx != b.midx)
                return a.midx < b.midx;
            return a.score < b.score;
        });

//    for (uint16_t sidx = 0; sidx < bds.get_num_seq(); ++ sidx)
//    {
//        for (auto const & motif : motifs)
//        {
//            uint16_t const hitidx = motif.uid * bds.get_num_seq() + sidx;
//            assert(hitidx < hits.size());
//            std::cerr << "Seq. " << sidx << ":\tFound " << hits[hitidx].size() << " hits for motif " << +motif.uid;
//            MotifScore sum = 0, min = 0, max = 0;
//            for (auto const & hit : hits[hitidx])
//            {
//                sum += std::get<2>(hit);
//                min = std::min(min, std::get<2>(hit));
//                max = std::max(max, std::get<2>(hit));
//            }
//            std::cerr << ", \tscore average " << sum / hits[hitidx].size() << "";
//            std::cerr << ", \tmin " << min << "";
//            std::cerr << ", \tmax " << max << "\n";
//        }
//    }

//    constexpr auto score_fn = [] (float score, size_t dist)
//    {
//        return score - (score * dist * dist / 900);
//    };

    #pragma omp parallel for num_threads(4)
    for (uint16_t sidx = 0u; sidx < bds.get_num_seq(); ++sidx)
    {
        std::vector<Hit> const & hitvec = hits[sidx];
        if (hitvec.empty())
            continue;

        auto left_end = hitvec.cbegin();
        auto right_end = left_end;

        do
        {
            std::vector<bool> diversity(num_motifs, false);
            while (right_end != hitvec.cend() && right_end->pos <= left_end->pos + 30)
            {
                diversity[right_end->midx] = true;
                ++right_end;
            }
            uint8_t div = 0;
            for (bool b : diversity)
                div += b ? 1 : 0;
            if (right_end - left_end > 1 && div > 1)
            {
                #pragma omp critical (print)
                {
                    std::cerr << "Result " << sidx << ": ";
                    for (auto it = left_end; it != right_end; ++it)
                        std::cerr << "(" << +it->midx << "|" << it->pos << "," << it->score << ") ";
                    std::cerr << std::endl;
                }
            }
            ++left_end;
            right_end = left_end;
        } while (right_end != hitvec.cend());
    }
}

} // namespace mars
