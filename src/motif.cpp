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

StemElement & StemloopMotif::new_stem()
{
    elements.emplace_back<StemElement>({});
    return std::get<StemElement>(elements.back());
}

LoopElement & StemloopMotif::new_loop(bool is_5prime)
{
    elements.emplace_back<LoopElement>({});
    auto & elem = std::get<LoopElement>(elements.back());
    elem.is_5prime = is_5prime;
    return elem;
}

std::vector<StemloopMotif> detect_stemloops(std::vector<int> const & bpseq, std::vector<int> const & plevel)
{
    struct PkInfo
    {
        int level;
        bool closing;
        Coordinate previous;
    };
    std::vector<PkInfo> pk_infos{};
    std::vector<StemloopMotif> stemloops;
    unsigned char id_cnt{0u};

    auto contains_loop = [&stemloops] (Coordinate const & outer)
    {
        return std::any_of(stemloops.crbegin(), stemloops.crend(), [&outer] (StemloopMotif const & inner)
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
    if (settings.exterior) // add long external and multiloops
    {
        MotifLen pos = 0;
        size_t const len = stemloops.size();
        for (size_t idx = 0; idx < len; ++idx) // no iterators, because we modify the vector
        {
            if (stemloops[idx].bounds.first > pos + 19u)
                stemloops.emplace_back(id_cnt++, Coordinate{pos, stemloops[idx].bounds.first - 1u});
            pos = stemloops[idx].bounds.second + 1;
        }
        if (bpseq.size() > pos + 19u || stemloops.empty())
            stemloops.emplace_back(id_cnt, Coordinate{pos, bpseq.size() - 1u});
    }
    return stemloops;
}

std::vector<StemloopMotif> create_motifs()
{
    if (settings.alignment_file.empty())
        return {};
    else if (settings.alignment_file.extension().string().find("mmo") != std::string::npos)
        return restore_motifs(settings.alignment_file);

    // Read the alignment
    Msa msa = read_msa(settings.alignment_file);

    // Find the stem loops
    std::vector<StemloopMotif> motifs = detect_stemloops(msa.structure.first, msa.structure.second);

    // Create a structure motif for each stemloop
    std::vector<std::future<void>> futures;
    for (StemloopMotif & motif : motifs)
        futures.push_back(pool->submit(&StemloopMotif::analyze, &motif, msa));

    for (auto & future : futures)
        future.wait();

    logger(1, "Found " << motifs.size() << " stem loops <== " << settings.alignment_file << std::endl);
    for (auto const & motif : motifs)
        logger(2, motif << std::endl);
    return motifs;
}

// private helper functions for analyze()
void check_gaps(int & current_gap, std::vector<std::unordered_map<MotifLen, SeqNum>> & gaps, int col, bool is_gap)
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

void filter_gaps(std::vector<std::unordered_map<MotifLen, SeqNum>> & gaps, SeqNum depth)
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
    if (queue.size() < 2)
        return;
    auto ptr = queue.cbegin();
    while (ptr != queue.cend() && ptr->first < log2f(settings.prune/100.f))
        ++ptr;
    if (ptr == queue.cend()) // avoid empty profile
        --ptr;
    queue.erase(queue.cbegin(), ptr); // erase if ratio < p %
}

void StemloopMotif::analyze(Msa const & msa)
{
    depth = msa.sequences.size();
    std::valarray<MotifLen> motif_len_stat(static_cast<MotifLen>(0), depth);
    std::vector<int> const & bpseq = msa.structure.first;

    auto make_stem = [this, &msa, &bpseq, &motif_len_stat] (int & left, int & right)
    {
        StemElement & elem = new_stem();
        std::vector<int> gap_stat(depth, -1);
        std::valarray<MotifLen> len_stat(static_cast<MotifLen>(0), depth);
        do
        {
            assert(bpseq[right] == left);
            elem.gaps.emplace_back();
            profile_char<bi_alphabet<seqan3::rna4>> prof{};
            for (auto &&[current_gap, len, seq] : seqan3::views::zip(gap_stat, len_stat, msa.sequences))
            {
                bool is_gap = prof.increment(seq[left], seq[right]);
                if (seq[left] != seqan3::gap())
                    ++len;
                if (seq[right] != seqan3::gap())
                    ++len;
                check_gaps(current_gap, elem.gaps, elem.prio.size(), is_gap);
            }
            elem.prio.push_back(prof.priority(depth));
            filter_profile(elem.prio.back());
            ++left;
            --right;
        }
        while (bpseq[left] == right);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps, elem.prio.size(), false);

        filter_gaps(elem.gaps, depth);
        elem.length = {len_stat.min(), len_stat.max(), len_stat.sum() / static_cast<float>(depth)};
        motif_len_stat += len_stat;
    };

    auto make_loop = [this, &msa, &bpseq, &motif_len_stat] (int & bpidx, bool is_5prime)
    {
        LoopElement & elem = new_loop(is_5prime);
        std::vector<int> gap_stat(depth, -1);
        std::valarray<MotifLen> len_stat(static_cast<MotifLen>(0), depth);
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
            elem.prio.push_back(prof.priority(depth));
            filter_profile(elem.prio.back());
            bpidx += (elem.is_5prime ? 1 : -1);
        }
        while ((bpseq[bpidx] < bounds.first || bpseq[bpidx] > bounds.second) &&
               bpidx >= bounds.first && bpidx <= bounds.second);

        for (int current_gap : gap_stat)
            check_gaps(current_gap, elem.gaps, elem.prio.size(), false);

        filter_gaps(elem.gaps, depth);
        elem.length = {len_stat.min(), len_stat.max(), len_stat.sum() / static_cast<float>(depth)};
        motif_len_stat += len_stat;
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
    length = {motif_len_stat.min(), motif_len_stat.max(), motif_len_stat.sum() / static_cast<float>(depth)};
}

void StemloopMotif::print_rssp(std::ofstream & os) const
{
    os << ">RSSP" << +uid << "|startpos=" << bounds.first << "|weight=1\n";
    std::deque<char> sequence{};
    std::deque<char> structure{};

    for (auto elemIt = elements.crbegin(); elemIt != elements.crend(); ++elemIt)
    {
        std::visit([&sequence, &structure] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), LoopElement>)
            {
                auto push = [&sequence, &structure, &element] (char c)
                {
                    if (element.is_5prime)
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
                for (auto profIt = element.prio.crbegin(); profIt != element.prio.crend(); ++profIt)
                {
                    if (profIt->size() > 1)
                        push('N');
                    else if (profIt->size() == 1)
                        push(profIt->cbegin()->second.to_char());
                }
            }
            else // stem
            {
                for (auto profIt = element.prio.crbegin(); profIt != element.prio.crend(); ++profIt)
                {
                    if (profIt->size() > 1)
                    {
                        sequence.push_front('N');
                        sequence.push_back('N');
                        structure.push_front('(');
                        structure.push_back(')');
                    }
                    else if (profIt->size() == 1)
                    {
                        std::pair<char, char> chrs = profIt->cbegin()->second.to_chars();
                        sequence.push_front(chrs.first);
                        sequence.push_back(chrs.second);
                        structure.push_front('(');
                        structure.push_back(')');
                    }
                }
            }
        }, *elemIt);
    }

    for (char const c : sequence)
        os << c;
    os << '\n';
    for (char const c : structure)
        os << c;
    os << '\n';
}

std::ostream & operator<<(std::ostream & os, StemloopMotif const & motif)
{
    os << "[" << +motif.uid << "] MOTIF pos = (" << motif.bounds.first << ".."
       << motif.bounds.second << "), len = (" << motif.length.min << ".." << motif.length.max << "), avg= "
       << std::fixed << std::setprecision(1) << motif.length.mean << "\n";
    for (auto elem = motif.elements.crbegin(); elem != motif.elements.crend(); ++elem)
    {
        std::visit([&os] (auto element)
        {
            if constexpr (std::is_same_v<decltype(element), LoopElement>)
                os << "\tLoop " << (element.is_5prime ? "5' " : "3' ");
            else
                os << "\tStem ";

            os << "L=(" << element.length.min << ".." << element.length.max << "), avg= "
               << std::fixed << std::setprecision(1) << element.length.mean << ":\t";

            for (int idx = element.prio.size() - 1; idx >= 0 ; --idx)
            {
                if (!element.prio[idx].empty())
                {
                    auto iter = element.prio[idx].crbegin();
                    if constexpr (std::is_same_v<decltype(element), StemElement>)
                    {
                        auto [ch1, ch2] = iter->second.to_chars();
                        os << "(" << ch1 << ch2;
                        for (++iter; iter != element.prio[idx].crend(); ++iter)
                        {
                            auto [ch3, ch4] = iter->second.to_chars();
                            os << "," << ch3 << ch4;
                        }
                    }
                    else
                    {
                        os << "(" << iter->second.to_char();
                        for (++iter; iter != element.prio[idx].crend(); ++iter)
                            os << "," << iter->second.to_char();
                    }
                    if (!element.gaps[idx].empty())
                        os << ",";
                }
                else
                    os << "(";
                for (auto const & gap : element.gaps[idx])
                    os << gap.first << "-";
                os << ") ";
            }
            os << "\n";
        }, *elem);
     }
    return os;
}

void store_rssp(std::vector<StemloopMotif> const & motifs)
{
    if (settings.structator_file.empty() || motifs.empty())
        return;

    std::ofstream ofs(settings.structator_file);
    if (ofs)
    {
        for (auto const & motif : motifs)
            motif.print_rssp(ofs);
        logger(1, "Stored " << motifs.size() << " motifs ==> " << settings.structator_file << std::endl);
    }
    ofs.close();
}

#if SEQAN3_WITH_CEREAL
std::vector<StemloopMotif> restore_motifs(std::filesystem::path const & motif_file)
{
    std::vector<StemloopMotif> motifs{};
    std::ifstream ifs{motif_file, std::ios::binary};
    if (ifs.good())
    {
        cereal::BinaryInputArchive iarchive{ifs};
        std::string version;
        iarchive(version);
        if (version[0] == '1')
        {
            iarchive(motifs);
            logger(1, "Restored " << motifs.size() << " motifs <== " << motif_file << std::endl);
        }
    }
    ifs.close();
    return motifs;
}

void store_motifs(std::vector<StemloopMotif> const & motifs)
{
    if (settings.motif_file.empty() || motifs.empty())
        return;
    std::ofstream ofs{settings.motif_file, std::ios::binary};
    if (ofs)
    {
        // Write the index to disk, including a version string.
        cereal::BinaryOutputArchive oarchive{ofs};
        std::string const version{"1 mars vector<StemloopMotif>\n"};
        oarchive(version);
        oarchive(motifs);
        logger(1, "Stored " << motifs.size() << " motifs ==> " << settings.motif_file << std::endl);
    }
    ofs.close();
}
#endif

} // namespace mars
