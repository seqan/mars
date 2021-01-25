#include <seqan3/std/ranges>

#include <seqan3/range/views/deep.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/zip.hpp>

#include "bi_alphabet.hpp"
#include "motif.hpp"

namespace mars
{

stemloop_type detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel)
{
    struct pk_info
    {
        int level;
        bool closing;
        std::pair<int, int> previous;
    };
    std::vector<pk_info> pk_infos{};
    stemloop_type stemloops;

    // 0-based indices
    for (auto &&[idx, bp, pk] : seqan3::views::zip(std::ranges::views::iota(0), bpseq, plevel))
    {
        if (pk == -1) // skip unpaired
            continue;

        while (pk + 1 > pk_infos.size()) // allocate a new pseudoknot layer
            pk_infos.push_back({0, false, {0, 0}});

        pk_info & status = pk_infos[pk];
        if (bp < idx) // close an interaction
        {
            status.previous = {bp, idx};
            status.closing = true;
            if (--status.level == 0)
                stemloops.push_back(status.previous);
        }
        else if (status.closing) // open an interaction (after closing the previous)
        {
            if (status.level > 0)
                stemloops.push_back(status.previous);
            status.level = 1;
            status.closing = false;
        }
        else // open another interaction
        {
            ++status.level;
        }
    }
    return std::move(stemloops);
}

stem_loop_motif analyze_stem_loop(msa_type const & msa, std::vector<int> const & bpseq, std::pair<int, int> const & pos)
{
    stem_loop_motif motif{};
    motif.bounds = pos;
    int left = pos.first;
    int right = pos.second;

    while (left <= right)
    {
        if (bpseq[left] == right) // stem
        {
            std::vector<int> gap_stat(msa.sequences.size(), -1);
            stem_element & stem = motif.new_stem();
            uint16_t col = 0u;
            do
            {
                assert(bpseq[right] == left);
                stem.gaps.emplace_back();
                profile_char<mars::bi_alphabet<seqan3::rna4>> prof{};
                for (auto &&[current_gap, seq] : seqan3::views::zip(gap_stat, msa.sequences))
                {
                    bool is_gap = prof.increment(seq[left], seq[right]);
                    if (is_gap && current_gap == -1)
                    {
                        current_gap = col;
                    }
                    else if (!is_gap && current_gap > -1)
                    {
                        auto [iter, succ] = stem.gaps[current_gap].emplace(col-current_gap, 1);
                        if (!succ)
                            ++(iter->second);
                        current_gap = -1;
                    }
                }
                stem.profile.push_back(prof);
                ++col;
                ++left;
                --right;
            }
            while (bpseq[left] == right);

            for (int current_gap : gap_stat)
            {
                if (current_gap > -1)
                {
                    auto[iter, succ] = stem.gaps[current_gap].emplace(col - current_gap, 1);
                    if (!succ)
                        ++(iter->second);
                }
            }
        }

        if (bpseq[left] < pos.first || bpseq[left] > pos.second) // 5' loop
        {
            std::vector<int> gap_stat(msa.sequences.size(), -1);
            loop_element & loop = motif.new_loop();
            loop.is_5prime = true;
            uint16_t col = 0u;
            do
            {
                loop.gaps.emplace_back();
                profile_char<seqan3::rna4> prof{};
                for (auto &&[current_gap, seq] : seqan3::views::zip(gap_stat, msa.sequences))
                {
                    bool is_gap = prof.increment(seq[left]);
                    if (is_gap && current_gap == -1)
                    {
                        current_gap = col;
                    }
                    else if (!is_gap && current_gap > -1)
                    {
                        auto [iter, succ] = loop.gaps[current_gap].emplace(col-current_gap, 1);
                        if (!succ)
                            ++(iter->second);
                        current_gap = -1;
                    }
                }
                loop.profile.push_back(prof);
                ++left;
                ++col;
            }
            while (bpseq[left] < pos.first || bpseq[left] > pos.second);

            for (int current_gap : gap_stat)
            {
                if (current_gap > -1)
                {
                    auto[iter, succ] = loop.gaps[current_gap].emplace(col - current_gap, 1);
                    if (!succ)
                        ++(iter->second);
                }
            }
        }
        else if (bpseq[right] < pos.first || bpseq[right] > pos.second) // 3' loop
        {
            std::vector<int> gap_stat(msa.sequences.size(), -1);
            loop_element & loop = motif.new_loop();
            loop.is_5prime = false;
            uint16_t col = 0u;
            do
            {
                loop.gaps.emplace_back();
                profile_char<seqan3::rna4> prof{};
                for (auto &&[current_gap, seq] : seqan3::views::zip(gap_stat, msa.sequences))
                {
                    bool is_gap = prof.increment(seq[right]);
                    if (is_gap && current_gap == -1)
                    {
                        current_gap = col;
                    }
                    else if (!is_gap && current_gap > -1)
                    {
                        auto [iter, succ] = loop.gaps[current_gap].emplace(col-current_gap, 1);
                        if (!succ)
                            ++(iter->second);
                        current_gap = -1;
                    }
                }
                loop.profile.push_back(prof);
                --right;
                ++col;
            }
            while (bpseq[right] < pos.first || bpseq[right] > pos.second);

            for (int current_gap : gap_stat)
            {
                if (current_gap > -1)
                {
                    auto[iter, succ] = loop.gaps[current_gap].emplace(col - current_gap, 1);
                    if (!succ)
                        ++(iter->second);
                }
            }
        }
        else
        {
            assert(false);
        }
    }

    return std::move(motif);
}

} // namespace mars
