#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/search.hpp>


#include "index.hpp"

namespace mars
{

void BiDirectionalIndex::create(std::filesystem::path const & filepath)
{
    if (filepath.empty())
        return;

    std::filesystem::path indexpath = filepath;
    indexpath += ".marsindex";

    // Check whether an index already exists.
    if (read_index(index, indexpath))
        return;

    // No index found: read genome and create an index.
    if (std::filesystem::exists(filepath))
    {
        // Generate the BiFM index.
        auto seqs = read_genome(filepath);
        index = Index{seqs};
        write_index(index, indexpath);
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

void BiDirectionalIndex::append_loop(std::pair<float, seqan3::rna4> item, bool left)
{
    assert(!queries.empty());
    seqan3::dna4_vector elem = queries.back();
    if (left)
        elem.insert(elem.begin(), item.second);
    else
        elem.push_back(item.second);
    queries.push_back(elem);
    scores.push_back(scores.back() + item.first);
//    seqan3::debug_stream << "add loop " << scores << "\t" << queries.back() << "\t" << item.first << "\n";
}

void BiDirectionalIndex::append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item)
{
    using seqan3::get;
    assert(!queries.empty());
    seqan3::dna4_vector elem = queries.back();
    seqan3::rna4 c = get<0>(stem_item.second);
    elem.insert(elem.begin(), c);
    c = get<1>(stem_item.second);
    elem.push_back(c);
    queries.push_back(elem);
    scores.push_back(scores.back() + stem_item.first);
//    seqan3::debug_stream << "add stem " << scores << "\t" << queries.back() << "\t" << stem_item.first << "\n";
}

void BiDirectionalIndex::backtrack()
{
    assert(!queries.empty());

    queries.pop_back();
    scores.pop_back();
}

float BiDirectionalIndex::get_score() const
{
    return scores.back();
}

bool BiDirectionalIndex::xdrop() const
{
    if (scores.size() < xdrop_dist)
        return false;
    else
        return scores.back() < scores[scores.size() - xdrop_dist];
}

size_t BiDirectionalIndex::compute_matches()
{
    assert(!queries.empty());
    assert(!index.empty());
    auto res = seqan3::search(queries.back(), index);
    matches.clear();
    for (auto & result : res)
        matches.emplace_back(result.reference_id(), result.reference_begin_position());
    return matches.size();
}

} // namespace mars
