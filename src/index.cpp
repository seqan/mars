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
    if (read_index(index, names, indexpath))
    {
        cursors.emplace_back(index);
        if (verbose > 0)
            std::cerr << "Using existing index <== " << indexpath << std::endl;
        return;
    }

    // No index found: read genome and create an index.
    if (std::filesystem::exists(filepath))
    {
        if (verbose > 0)
            std::cerr << "Read genome <== " << filepath << std::endl;
        std::vector<seqan3::dna4_vector> seqs{};
        read_genome(seqs, names, filepath);
        // Generate the BiFM index.
        index = Index{seqs};
        cursors.emplace_back(index);
        write_index(index, names, indexpath);
        if (verbose > 0)
            std::cerr << "Created index ==> " << indexpath << std::endl;
    }
    else
    {
        std::ostringstream err_msg{};
        err_msg << "Could not find the genome file <== " << filepath << "[.marsindex";
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
    {
        assert(seq < hits.size());
        hits[seq].emplace_back(pos + max_offset - motif.bounds.first,
                               motif.uid,
                               scores.back());
    }
}

} // namespace mars
