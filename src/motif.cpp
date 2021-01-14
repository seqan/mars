#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
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
    auto ali = msa.sequences | seqan3::views::deep{seqan3::views::slice}(pos.first, pos.second + 1);
    auto pairs = bpseq | seqan3::views::slice(pos.first, pos.second + 1);
    int left = pos.first;
    int right = pos.second;

    while (left <= right)
    {
        if (bpseq[left] == right) // stem
        {
            stem_element & stem = motif.new_stem();
            do
            {
                assert(bpseq[right] == left);
                seqan3::debug_stream << "add to stack: " << left << " & " << right;
                profile_char<mars::bi_alphabet<seqan3::rna4>> p{};
                for (auto const & seq : msa.sequences)
                {
                    bool is_gap = p.increment(seq[left], seq[right]);
                }
                seqan3::debug_stream << "\t profile " << p << std::endl;
                stem.profile.push_back(p);
                ++left;
                --right;
            }
            while (bpseq[left] == right);
        }

        if (bpseq[left] < pos.first || bpseq[left] > pos.second) // 5' loop
        {
            loop_element & loop = motif.new_loop();
            loop.is_5prime = true;
            do
            {
                seqan3::debug_stream << "left loop: " << left;
                profile_char<seqan3::rna4> p{};
                for (auto const & seq : msa.sequences)
                {
                    bool is_gap = p.increment(seq[left]);
                }
                seqan3::debug_stream << "\t profile " << p << std::endl;
                loop.profile.push_back(p);
                ++left;
            }
            while (bpseq[left] < pos.first || bpseq[left] > pos.second);
        }
        else if (bpseq[right] < pos.first || bpseq[right] > pos.second) // 3' loop
        {
            loop_element & loop = motif.new_loop();
            loop.is_5prime = false;
            do
            {
                seqan3::debug_stream << "right loop: " << right;
                profile_char<seqan3::rna4> p{};
                for (auto const & seq : msa.sequences)
                {
                    bool is_gap = p.increment(seq[right]);
                }
                seqan3::debug_stream << "\t profile " << p << std::endl;
                loop.profile.push_back(p);
                --right;
            }
            while (bpseq[right] < pos.first || bpseq[right] > pos.second);
        }
        else
        {
            assert(false);
        }
    }

    return std::move(motif);
}

} // namespace mars
