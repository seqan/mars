#include <seqan3/range/views/zip.hpp>

#ifdef MARS_WITH_OPENMP
    #include <omp.h>
#endif

#include "search.hpp"
#include "settings.hpp"

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
    if (mars::verbose > 0)
        std::cerr << "Start the motif search...";
    assert(motifs.size() <= UINT8_MAX);
    uint8_t const num_motifs = motifs.size();
    hits.resize(bds.number_of_seq());

    for (auto const & motif : motifs)
        bds.update_max_offset(motif.bounds.first);

//    #pragma omp parallel for num_threads(2)
    for (uint8_t midx = 0; midx < num_motifs; ++midx)
    {
        auto const & motif = motifs[midx];
        // start with the hairpin
        auto const iter = motif.elements.crbegin();
        recurse_search<LoopElement>(motif, iter, 0);
        if (verbose > 0)
            std::cerr << "  " << (100*(midx+1)/num_motifs) << "%";
    }
    if (verbose > 0)
        std::cerr << std::endl;

    #pragma omp parallel for num_threads(4)
    for (size_t sidx = 0u; sidx < bds.number_of_seq(); ++sidx)
        std::sort(hits[sidx].begin(), hits[sidx].end(), [] (Hit const & a, Hit const & b)
        {
            if (a.pos != b.pos)
                return a.pos < b.pos;
            if (a.midx != b.midx)
                return a.midx < b.midx;
            return a.score < b.score;
        });

    size_t num_results = 0;
//    #pragma omp parallel for num_threads(4)
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
                hit.score = static_cast<float>(std::max(0.0, hit.score - (0.05 * pos_diff * pos_diff)));
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

            if (diversity > 1)
            {
                locations.emplace(hit_score, diversity, base_pos, sidx);
//                std::cout << ">" << std::left << std::setw(35) << bds.get_name(sidx) << "\t" << sidx << "\t"
//                          << base_pos << "\t" << +diversity << "\t" << hit_score << std::endl;
                ++num_results;
            }

            left_end = right_end;
        } while (right_end != hitvec.cend());
    }
    if (verbose > 1)
        std::cerr << "Found " << num_results << " matches." << std::endl;
}

} // namespace mars
