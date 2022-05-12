// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <seqan3/std/ranges>
#include <set>
#include <tuple>
#include <valarray>

#include <seqan3/utility/views/deep.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#endif

#include "motif.hpp"
#include "settings.hpp"

namespace mars
{

// private helper functions for analyze()
void check_gaps(int & current_gap, std::vector<std::unordered_map<Position, size_t>> & gaps, int col, bool is_gap)
{
    if (is_gap && current_gap == -1) // gap start
    {
        current_gap = col;
    }
    else if (!is_gap && current_gap > -1) // gap end
    {
        auto[iter, succ] = gaps[col - 1].emplace(col - current_gap, 1);
        if (!succ)
            ++(iter->second);
        current_gap = -1;
    }
}

void filter_gaps(std::vector<std::unordered_map<Position, size_t>> & gaps, size_t depth)
{
    for (auto & map : gaps)
    {
        auto iter = map.begin();
        while (iter != map.end())
        {
            if (iter->second * 200 <= depth * settings.prune) // erase if gap ratio < p/2 %
                iter = map.erase(iter);
            else
                ++iter;
        }
    }
}

void filter_profile(auto & queue)
{
    if (queue.size() < 2 || settings.prune == 0)
        return;
    auto ptr = queue.cbegin();
    while (ptr != queue.cend() && ptr->first < log2f(settings.prune/100.f))
        ++ptr;
    if (ptr == queue.cend()) // avoid empty profile
        --ptr;
    queue.erase(queue.cbegin(), ptr); // erase if ratio < p %
}

void Stemloop::analyze(Msa const & msa)
{
    depth = msa.sequences.size();
    std::valarray<Position> sl_len_stat(static_cast<Position>(0), depth);
    std::vector<int> const & bpseq = msa.structure.first;

    auto make_stem = [this, &msa, &bpseq, &sl_len_stat] (int & left, int & right)
    {
        StemElement & elem = std::get<StemElement>(elements.emplace_back<StemElement>({}));
        std::vector<int> gap_stat(depth, -1);
        std::valarray<Position> len_stat(static_cast<Position>(0), depth);
        do
        {
            assert(bpseq[right] == left);
            elem.gaps.emplace_back();
            profile_char<bi_alphabet<seqan3::gapped<seqan3::rna4>>> prof{};
            for (auto &&[current_gap, len, seq] : seqan3::views::zip(gap_stat, len_stat, msa.sequences))
            {
                bool is_gap = seq[left] == seqan3::gap() && seq[right] == seqan3::gap();
                if (!is_gap)
                {
                    prof.increment(seq[left], seq[right]);
                    len += (seq[left] != seqan3::gap() && seq[right] != seqan3::gap()) ? 2 : 1;
                }
                check_gaps(current_gap, elem.gaps, elem.prio.size(), is_gap);
            }
            elem.prio.push_back(priority(prof, depth));
            filter_profile(elem.prio.back());
            ++left;
            --right;
        }
        while (bpseq[left] == right);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps, elem.prio.size(), false);

        filter_gaps(elem.gaps, depth);
        sl_len_stat += len_stat;
        std::reverse(elem.gaps.begin(), elem.gaps.end());
        std::reverse(elem.prio.begin(), elem.prio.end());
    };

    auto make_loop = [this, &msa, &bpseq, &sl_len_stat] (int & bpidx, bool leftsided)
    {
        LoopElement & elem = std::get<LoopElement>(elements.emplace_back<LoopElement>({}));
        elem.leftsided = leftsided;
        std::vector<int> gap_stat(depth, -1);
        std::valarray<Position> len_stat(static_cast<Position>(0), depth);
        do
        {
            elem.gaps.emplace_back();
            profile_char<seqan3::rna4> prof{};
            for (auto &&[current_gap, len, seq] : seqan3::views::zip(gap_stat, len_stat, msa.sequences))
            {
                bool is_gap = prof.increment(seq[bpidx]);
                if (!is_gap)
                    ++len;
                check_gaps(current_gap, elem.gaps, elem.prio.size(), is_gap);
            }
            elem.prio.push_back(priority(prof, depth));
            filter_profile(elem.prio.back());
            bpidx += (elem.leftsided ? 1 : -1);
        }
        while ((bpseq[bpidx] < bounds.first || bpseq[bpidx] > bounds.second) &&
               bpidx >= bounds.first && bpidx <= bounds.second);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps, elem.prio.size(), false);

        filter_gaps(elem.gaps, depth);
        sl_len_stat += len_stat;
        std::reverse(elem.gaps.begin(), elem.gaps.end());
        std::reverse(elem.prio.begin(), elem.prio.end());
    };

    int left = bounds.first;
    int right = bounds.second;
    while (left <= right)
    {
        if (bpseq[left] == right) // stem
            make_stem(left, right);
        else if (bpseq[right] < bounds.first || bpseq[right] > bounds.second) // 3' loop
            make_loop(right, false);
        else if (bpseq[left] < bounds.first || bpseq[left] > bounds.second) // 5' loop
            make_loop(left, true);
        else
        {
            logger(0, "Unexpected condition! " << bpseq << std::endl);
            throw std::runtime_error("The structure is inconsistent."); // prevent endless loop
        }
    }
    length = {sl_len_stat.min(), sl_len_stat.max()};
    std::reverse(elements.begin(), elements.end());
}

void Stemloop::print_rssp(std::ofstream & os) const
{
    os << ">RSSP" << +uid << "|startpos=" << bounds.first << "|weight=1\n";
    std::deque<char> sequence{};
    std::deque<char> structure{};

    for (auto const & elem : elements)
    {
        std::visit([&sequence, &structure] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), LoopElement>)
            {
                auto push = [&sequence, &structure, &element] (char c)
                {
                    if (element.leftsided)
                    {
                        sequence.push_front(c);
                        structure.push_front('.');
                    }
                    else
                    {
                        sequence.push_back(c);
                        structure.push_back('.');
                    }
                };
                for (auto const & prof : element.prio)
                {
                    if (prof.size() > 1)
                        push('N');
                    else if (prof.size() == 1)
                        push(prof.front().second.to_char());
                }
            }
            else // stem
            {
                for (auto const & prof : element.prio)
                {
                    if (prof.size() > 1)
                    {
                        sequence.push_front('N');
                        sequence.push_back('N');
                        structure.push_front('(');
                        structure.push_back(')');
                    }
                    else if (prof.size() == 1)
                    {
                        std::pair<char, char> chrs = prof.front().second.to_chars();
                        sequence.push_front(chrs.first);
                        sequence.push_back(chrs.second);
                        structure.push_front('(');
                        structure.push_back(')');
                    }
                }
            }
        }, elem);
    }

    for (char const c : sequence)
        os << c;
    os << '\n';
    for (char const c : structure)
        os << c;
    os << '\n';
}

std::ostream & operator<<(std::ostream & os, Stemloop const & stemloop)
{
    os << "[" << (+stemloop.uid + 1) << "] STEMLOOP pos = (" << stemloop.bounds.first << ".."
       << stemloop.bounds.second << "), len = (" << stemloop.length.first << ".." << stemloop.length.second << ")\n";
    for (auto const & elem : stemloop.elements)
    {
        std::visit([&os] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), LoopElement>)
                os << "\tLoop " << (element.leftsided ? "5' " : "3' ");
            else
                os << "\tStem    ";

            for (auto && [prio, gaps] : seqan3::views::zip(element.prio, element.gaps))
            {
                if (!prio.empty())
                {
                    auto iter = prio.crbegin();
                    if constexpr (std::is_same_v<decltype(element), StemElement>)
                    {
                        auto [ch1, ch2] = iter->second.to_chars();
                        os << "(" << ch1 << ch2;
                        for (++iter; iter != prio.crend(); ++iter)
                        {
                            auto [ch3, ch4] = iter->second.to_chars();
                            os << "," << ch3 << ch4;
                        }
                    }
                    else
                    {
                        os << "(" << iter->second.to_char();
                        for (++iter; iter != prio.crend(); ++iter)
                            os << "," << iter->second.to_char();
                    }
                    if (!gaps.empty())
                        os << ",";
                }
                else
                    os << "(";
                for (auto const & gap : gaps)
                    os << gap.first << "-";
                os << ") ";
            }
            os << "\n";
        }, elem);
    }
    return os;
}

void store_rssp(Motif const & motif)
{
    if (settings.structator_file.empty() || motif.empty())
        return;

    std::ofstream ofs(settings.structator_file);
    if (ofs)
    {
        for (auto const & stemloop : motif)
            stemloop.print_rssp(ofs);
        logger(1, "Stored " << motif.size() << " stemloops ==> " << settings.structator_file << std::endl);
    }
    ofs.close();
}

Motif create_motif()
{
    if (settings.alignment_file.empty())
        return {};
#if SEQAN3_WITH_CEREAL
    else if (settings.alignment_file.extension().string().find("mmo") != std::string::npos)
        return restore_motif(settings.alignment_file);
#endif

    // Read the alignment
    Msa msa = read_msa(settings.alignment_file);

    // Find the stem loops
    Motif motif = detect_stemloops(msa.structure.first, msa.structure.second);

    // Analyze each stemloop
    std::vector<std::future<void>> futures;
    for (Stemloop & stemloop : motif)
        futures.push_back(pool->submit(&Stemloop::analyze, &stemloop, msa));

    for (auto & future : futures)
        future.wait();

    logger(1, "Found " << motif.size() << " stemloops <== " << settings.alignment_file << std::endl);
    for (auto const & stemloop : motif)
    logger(2, stemloop << std::endl);
    return motif;
}

Motif detect_stemloops(std::vector<int> const & bpseq, std::vector<int> const & plevel)
{
    struct PkInfo
    {
        int level;
        bool closing;
        Bounds previous;
    };
    std::vector<PkInfo> pk_infos{};
    Motif stemloops;
    unsigned char id_cnt{0u};

    auto contains_loop = [&stemloops] (Bounds const & outer)
    {
        return std::any_of(stemloops.crbegin(), stemloops.crend(), [&outer] (Stemloop const & inner)
        {
            return outer.first < inner.bounds.first && inner.bounds.second < outer.second;
        });
    };

    // 0-based indices
    for (auto &&[idx, bp, pk] : seqan3::views::zip(std::ranges::views::iota(0), bpseq, plevel))
    {
        if (pk == -1) // skip unpaired
            continue;

        while (pk + 1 > static_cast<int>(pk_infos.size())) // allocate a new pseudoknot layer
            pk_infos.push_back({0, false, {0, 0}});

        PkInfo & status = pk_infos[pk];
        if (bp < idx) // close an interaction
        {
            status.previous = {bp, idx};
            status.closing = true;
            if (--status.level == 0 && !contains_loop(status.previous))
                stemloops.emplace_back(id_cnt++, status.previous);
        }
        else if (status.closing) // open an interaction (after closing the previous)
        {
            if (status.level > 0 && !contains_loop(status.previous))
                stemloops.emplace_back(id_cnt++, status.previous);
            status.level = 1;
            status.closing = false;
        }
        else // open another interaction
        {
            ++status.level;
        }
    }
    if (!settings.limit) // add long external and multiloops
    {
        Position pos = 0;
        size_t const len = stemloops.size();
        for (size_t idx = 0; idx < len; ++idx) // no iterators, because we modify the vector
        {
            if (stemloops[idx].bounds.first > pos + 19u)
                stemloops.emplace_back(id_cnt++, Bounds{pos, stemloops[idx].bounds.first - 1u});
            pos = stemloops[idx].bounds.second + 1;
        }
        if (bpseq.size() > pos + 19u || stemloops.empty())
            stemloops.emplace_back(id_cnt, Bounds{pos, bpseq.size() - 1u});
    }
    return stemloops;
}

#if SEQAN3_WITH_CEREAL
Motif restore_motif(std::filesystem::path const & motif_file)
{
    Motif motif{};
    std::ifstream ifs{motif_file, std::ios::binary};
    if (ifs.good())
    {
        cereal::BinaryInputArchive iarchive{ifs};
        std::string version;
        iarchive(version);
        if (version[0] == '1')
        {
            iarchive(motif);
            logger(1, "Restored " << motif.size() << " stemloops <== " << motif_file << std::endl);
        }
    }
    ifs.close();
    return motif;
}

void store_motif(Motif const & motif)
{
    if (settings.motif_file.empty() || motif.empty())
        return;
    std::ofstream ofs{settings.motif_file, std::ios::binary};
    if (ofs)
    {
        // Write the index to disk, including a version string.
        cereal::BinaryOutputArchive oarchive{ofs};
        std::string const version{"1 mars vector<Stemloop>\n"};
        oarchive(version);
        oarchive(motif);
        logger(1, "Stored " << motif.size() << " stemloops ==> " << settings.motif_file << std::endl);
    }
    ofs.close();
}
#endif

} // namespace mars
