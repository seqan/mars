#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/search.hpp>

#include "index.hpp"
#include "settings.hpp"

namespace mars
{

void BiDirectionalIndex::create(std::filesystem::path const & filepath)
{
    if (filepath.empty())
        return;

    std::filesystem::path indexpath = filepath;
    indexpath += ".marsindex";

    // Check whether an index already exists.
    if (read_index(index, index_num_seq, indexpath))
    {
        cursors.emplace_back(index);
        if (verbose > 0)
            std::cerr << "Using existing index file: " << indexpath << std::endl;
        return;
    }

    // No index found: read genome and create an index.
    if (std::filesystem::exists(filepath))
    {
        // Generate the BiFM index.
        std::vector<std::vector<seqan3::dna4>> seqs = read_genome(filepath);
        index = Index{seqs};
        index_num_seq = seqs.size();
        cursors.emplace_back(index);
        write_index(index, index_num_seq, indexpath);
        if (verbose > 0)
            std::cerr << "Read genome from " << filepath << "\nCreate index file " << indexpath << std::endl;
    }
    else
    {
        std::ostringstream err_msg{};
        err_msg << "Could not find the genome file: " << filepath << "[.marsindex";
#ifdef SEQAN3_HAS_ZLIB
        err_msg << ".gz";
#endif
        err_msg << "]";
        throw seqan3::file_open_error(err_msg.str());
    }
}

bool BiDirectionalIndex::append_loop(std::pair<float, seqan3::rna4> item, bool left)
{
    bool succ;
    seqan3::bi_fm_index_cursor<Index> new_cur(cursors.back());

    if (left)
        succ = new_cur.extend_left(item.second);
    else
        succ = new_cur.extend_right(item.second);

    if (succ)
    {
        cursors.push_back(new_cur);
        scores.push_back(scores.back() + item.first);
    }

    return succ;
}

bool BiDirectionalIndex::append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item)
{
    seqan3::bi_fm_index_cursor<Index> new_cur(cursors.back());
    using seqan3::get;
    seqan3::rna4 c = get<0>(stem_item.second);
    bool succ = new_cur.extend_left(c);
    if (succ)
    {
        c = get<1>(stem_item.second);
        succ = new_cur.extend_right(c);
    }
    if (succ)
    {
        cursors.push_back(new_cur);
        scores.push_back(scores.back() + stem_item.first);
    }
    return succ;
}

void BiDirectionalIndex::backtrack()
{
    scores.pop_back();
    cursors.pop_back();
}

bool BiDirectionalIndex::xdrop() const
{
    if (scores.size() < xdrop_dist)
        return false;
    else
        return scores.back() < scores[scores.size() - xdrop_dist];
}

void BiDirectionalIndex::compute_hits(std::vector<std::vector<Hit>> & hits, StemloopMotif const & motif) const
{
    for (auto && [seq, pos] : cursors.back().locate())
        hits[seq].emplace_back(pos >= motif.bounds.first ? pos - motif.bounds.first : 0,
                               motif.uid,
                               scores.back());
}

} // namespace mars
