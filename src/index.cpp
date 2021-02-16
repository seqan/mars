#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/search.hpp>

#include "index.hpp"
#include "input_output.hpp"

namespace mars
{

Index create_index(std::filesystem::path const & filepath)
{
    if (filepath.empty())
        return {};

    std::filesystem::path indexpath = filepath;
    indexpath += ".marsindex";

    // Check whether an index already exists.
    if (std::filesystem::exists(indexpath))
    {
        std::ifstream ifs{indexpath, std::ios::binary};
        if (ifs.good())
        {
            cereal::BinaryInputArchive iarchive{ifs};

            // Verify the version string.
            std::string version;
            iarchive(version);
            assert(version[0] == '1');

            // Read the index from disk.
            Index index;
            iarchive(index);
            ifs.close();
            return std::move(index);
        }
        ifs.close();
    }

    // No index found: read genome and create an index.
    if (std::filesystem::exists(filepath))
    {
        // Generate the BiFM index.
        auto seqs = read_genome(filepath);
        Index index{seqs};
        std::ofstream ofs{indexpath, std::ios::binary};
        if (ofs)
        {
            // Write the index to disk, including a version string.
            cereal::BinaryOutputArchive oarchive{ofs};
            std::string version{"1 mars bi_fm_index rna5 collection\n"};
            oarchive(version);
            oarchive(index);
        }
        ofs.close();
        return std::move(index);
    }
    else
    {
        throw seqan3::file_open_error(std::string{"Could not find the genome file: "} + std::string{filepath});
    }
}

void BiDirectionalSearch::append_loop(std::pair<float, seqan3::rna4> item, bool left)
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

void BiDirectionalSearch::append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item)
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

void BiDirectionalSearch::backtrack()
{
    assert(!queries.empty());

    queries.pop_back();
    scores.pop_back();
}

float BiDirectionalSearch::get_score() const
{
    return scores.back();
}

bool BiDirectionalSearch::xdrop() const
{
    if (scores.size() < xdrop_dist)
        return false;
    else
        return scores.back() < scores[scores.size() - xdrop_dist];
}

size_t BiDirectionalSearch::compute_matches()
{
    assert(!queries.empty());
    auto res = seqan3::search(queries.back(), index);
    matches.clear();
    for (auto & result : res)
        matches.emplace_back(result.reference_id(), result.reference_begin_position());
    return matches.size();
}

} // namespace mars
