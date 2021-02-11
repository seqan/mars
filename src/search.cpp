#include <seqan3/range/views/zip.hpp>

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
void SearchGenerator::recurse_search(ElementIter const & it, MotifLen idx)
{
//    std::cerr << score << " == " << bds.get_score() << std::endl;
    if (bds.xdrop())
        return;

    auto const & elem = std::get<MotifElement>(*it);

    if (idx == elem.profile.size())
    {
        auto const next = it + 1;
        if (next != end)
        {
            if (std::holds_alternative<StemElement>(*next))
                recurse_search<StemElement>(next, 0);
            else
                recurse_search<LoopElement>(next, 0);
        }
        else
        {
            size_t num = bds.compute_matches();
            for (auto && [seq, pos] : bds.matches)
                hits.emplace_back(MotifNum{1}, seq, pos, bds.get_score());
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

        recurse_search<MotifElement>(it, idx + 1);
        bds.backtrack();
    }

    // try gaps
    for (auto && [len, num] : elem.gaps[idx])
        recurse_search<MotifElement>(it, idx + len);
}

void SearchGenerator::find_motif(StemloopMotif const & motif)
{
    std::cerr << "Let's search motif no. " << +motif.uid << "\n";

    // start with the hairpin
    auto const iter = motif.elements.crbegin();
    end = motif.elements.crend();
    recurse_search<LoopElement>(iter, 0);

    std::cerr << "Found " << hits.size() << " Matches";
    MotifScore sum = 0, min = 0, max = 0;
    for (auto const & hit : hits)
    {
        sum += hit.score;
        min = std::min(min, hit.score);
        max = std::max(max, hit.score);
    }
    std::cerr << ", score average " << sum/hits.size() << "";
    std::cerr << ", min " << min << "";
    std::cerr << ", max " << max << "\n";
}

} // namespace mars
