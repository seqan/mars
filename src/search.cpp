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
void SearchGenerator::recurse_search(MotifNum uid, ElementIter const & elem_it, MotifLen idx)
{
//    std::cerr << score << " == " << bds.get_score() << std::endl;
    if (bds.xdrop())
        return;

    auto const & elem = std::get<MotifElement>(*elem_it);

    if (idx == elem.profile.size())
    {
        auto const next = elem_it + 1;
        if (next != end_it[uid])
        {
            if (std::holds_alternative<StemElement>(*next))
                recurse_search<StemElement>(uid, next, 0);
            else
                recurse_search<LoopElement>(uid, next, 0);
        }
        else
        {
            size_t num = bds.compute_matches();
            for (auto && [seq, pos] : bds.matches)
                hits[uid].emplace_back(seq, pos, bds.get_score());
        }
        return;
    }

    auto prio = priority(elem.profile[idx]);

    // try to extend the pattern
    for (auto opt = prio.crbegin(); opt != prio.crend(); ++opt)
    {
        if constexpr (std::is_same_v<MotifElement, LoopElement>)
            bds.append_loop(*opt, false);
        else
            bds.append_stem(*opt);

        recurse_search<MotifElement>(uid, elem_it, idx + 1);
        bds.backtrack();
    }

    // try gaps
    for (auto && [len, num] : elem.gaps[idx])
        recurse_search<MotifElement>(uid, elem_it, idx + len);
}

void SearchGenerator::find_motifs(std::vector<StemloopMotif> const & motifs)
{
    hits.resize(motifs.size());
    end_it.resize(motifs.size());

    #pragma omp parallel for num_threads(2)
    for (size_t idx = 0; idx < motifs.size(); ++idx)
    {
        auto const & motif = motifs[idx];
        // start with the hairpin
        auto const iter = motif.elements.crbegin();
        end_it[motif.uid] = motif.elements.crend();
        recurse_search<LoopElement>(motif.uid, iter, 0);
    }

    for (auto const & motif : motifs)
    {
        std::cerr << "Found " << hits[motif.uid].size() << " matches for motif " << +motif.uid;
        MotifScore sum = 0, min = 0, max = 0;
        for (auto const & hit : hits[motif.uid])
        {
            sum += hit.score;
            min = std::min(min, hit.score);
            max = std::max(max, hit.score);
        }
        std::cerr << ", score average " << sum/hits[motif.uid].size() << "";
        std::cerr << ", min " << min << "";
        std::cerr << ", max " << max << "\n";
    }
}

} // namespace mars
