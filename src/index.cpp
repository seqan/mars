#include <seqan3/search/search.hpp>

#include "index.hpp"
#include "input_output.hpp"

namespace mars
{

index_type create_index(std::filesystem::path const & filepath)
{
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
            index_type index;
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
        index_type index{seqs};
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

void bi_directional_search::append_5prime(seqan3::dna4 left)
{
    assert(!queries.empty());
    seqan3::dna4_vector elem = queries.back();
    elem.insert(elem.begin(), left);
    queries.push_back(elem);
}

void bi_directional_search::append_3prime(seqan3::dna4 right)
{
    assert(!queries.empty());
    seqan3::dna4_vector elem = queries.back();
    elem.push_back(right);
    queries.push_back(elem);
}

void bi_directional_search::append_stem(seqan3::dna4 left, seqan3::dna4 right)
{
    assert(!queries.empty());
    seqan3::dna4_vector elem = queries.back();
    elem.insert(elem.begin(), left);
    elem.push_back(right);
    queries.push_back(elem);
}

bool bi_directional_search::backtrack()
{
    assert(!queries.empty());
    if (queries.back().empty())
        return false;

    queries.pop_back();
    return true;
}

size_t bi_directional_search::compute_matches()
{
    assert(!queries.empty());
    auto res = seqan3::search(queries.back(), index);
    matches.clear();
    for (auto & result : res)
        matches.emplace_back(result.reference_id(), result.reference_begin_position());
    return matches.size();
}

} // namespace mars
