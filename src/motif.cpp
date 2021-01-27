#include <seqan3/std/ranges>
#include <valarray>

#include <seqan3/range/views/deep.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/zip.hpp>

#include "motif.hpp"

namespace mars
{

stem_element & stemloop_motif::new_stem()
{
    elements.emplace_back<stem_element>({});
    return std::get<stem_element>(elements.back());
}

loop_element & stemloop_motif::new_loop(bool is_5prime)
{
    elements.emplace_back<loop_element>({});
    auto & elem = std::get<loop_element>(elements.back());
    elem.is_5prime = is_5prime;
    return elem;
}

std::vector<stemloop_type> detect_stem_loops(std::vector<int> const & bpseq, std::vector<int> const & plevel)
{
    struct pk_info
    {
        int level;
        bool closing;
        stemloop_type previous;
    };
    std::vector<pk_info> pk_infos{};
    std::vector<stemloop_type> stemloops;

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

// private helper function for analyze_stem_loop
void check_gaps(int & current_gap, std::unordered_map<uint16_t, uint16_t> & gaps, int col, bool is_gap)
{
    if (is_gap && current_gap == -1)
    {
        current_gap = col;
    }
    else if (!is_gap && current_gap > -1)
    {
        auto[iter, succ] = gaps.emplace(col - current_gap, 1);
        if (!succ)
            ++(iter->second);
        current_gap = -1;
    }
};

stemloop_motif analyze_stem_loop(msa_type const & msa, std::vector<int> const & bpseq, stemloop_type const & pos)
{
    stemloop_motif motif{};
    motif.bounds = pos;
    std::valarray<size_t> motif_len_stat(0ul, msa.sequences.size());

    auto make_stem = [&msa, &bpseq, &motif, &motif_len_stat] (int & left, int & right)
    {
        stem_element & elem = motif.new_stem();
        std::vector<int> gap_stat(msa.sequences.size(), -1);
        std::valarray<size_t> len_stat(0ul, msa.sequences.size());
        do
        {
            assert(bpseq[right] == left);
            elem.gaps.emplace_back();
            profile_char<mars::bi_alphabet<seqan3::rna4>> prof{};
            for (auto &&[current_gap, len, seq] : seqan3::views::zip(gap_stat, len_stat, msa.sequences))
            {
                bool is_gap = prof.increment(seq[left], seq[right]);
                if (seq[left] != seqan3::gap())
                    ++len;
                if (seq[right] != seqan3::gap())
                    ++len;
                check_gaps(current_gap, elem.gaps[current_gap], elem.profile.size(), is_gap);
            }
            elem.profile.push_back(prof);
            ++left;
            --right;
        }
        while (bpseq[left] == right);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps[current_gap], elem.profile.size(), false);

        elem.length = {len_stat.min(), len_stat.max(), len_stat.sum() * 1.f / msa.sequences.size()};
        motif_len_stat += len_stat;
    };

    auto make_loop = [&msa, &bpseq, &motif, &motif_len_stat] (int & bpidx, bool is_5prime)
    {
        loop_element & elem = motif.new_loop(is_5prime);
        std::vector<int> gap_stat(msa.sequences.size(), -1);
        std::valarray<size_t> len_stat(0ul, msa.sequences.size());
        do
        {
            elem.gaps.emplace_back();
            profile_char<seqan3::rna4> prof{};
            for (auto &&[current_gap, len, seq] : seqan3::views::zip(gap_stat, len_stat, msa.sequences))
            {
                bool is_gap = prof.increment(seq[bpidx]);
                if (!is_gap)
                    ++len;
                check_gaps(current_gap, elem.gaps[current_gap], elem.profile.size(), is_gap);
            }
            elem.profile.push_back(prof);
            bpidx += (elem.is_5prime ? 1 : -1);
        }
        while (bpseq[bpidx] < motif.bounds.first || bpseq[bpidx] > motif.bounds.second);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps[current_gap], elem.profile.size(), false);

        elem.length = {len_stat.min(), len_stat.max(), len_stat.sum() * 1.f / msa.sequences.size()};
        motif_len_stat += len_stat;
    };

    int left = motif.bounds.first;
    int right = motif.bounds.second;
    while (left <= right)
    {
        if (bpseq[left] == right) // stem
            make_stem(left, right);

        if (bpseq[left] < motif.bounds.first || bpseq[left] > motif.bounds.second) // 5' loop
            make_loop(left, true);
        else if (bpseq[right] < motif.bounds.first || bpseq[right] > motif.bounds.second) // 3' loop
            make_loop(right, false);
    }
    motif.length = {motif_len_stat.min(), motif_len_stat.max(), motif_len_stat.sum() * 1.f / msa.sequences.size()};
    return std::move(motif);
}

std::ostream & operator<<(std::ostream & os, stemloop_motif const & motif)
{
    os << "Motif pos = (" << motif.bounds.first << ", " << motif.bounds.second << "), "
       << "len = (" << motif.length.min << ", " << motif.length.max << ", " << motif.length.mean << ")\n";
    for (auto const & el : motif.elements)
    {
        std::visit([&os] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), mars::loop_element>)
                os << "Loop " << (element.is_5prime ? "5' " : "3' ");
            else
                os << "Stem ";

            os << "L=(" << element.length.min << ", " << element.length.max << ", " << element.length.mean << ") :\t";

            for (auto const & profile_char : element.profile)
                os << profile_char << ' ';
            os << "\nGaps: ";
            for (int idx = 0; idx < element.gaps.size(); ++idx)
                if (!element.gaps[idx].empty())
                {
                    os << "\t" << idx << ": ";
                    for (auto && [key, val] : element.gaps[idx])
                        os << "(" << key << "," << val << ")";
                }
            os << "\n";
        }, el);
     }
    return os;
}

} // namespace mars
