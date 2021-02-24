#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/search.hpp>

#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/gz_istream.hpp>
    #include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

#include "index.hpp"
#include "input_output.hpp"

namespace mars
{

void BiDirectionalSearch::write_index(std::filesystem::path & indexpath)
{
#ifdef SEQAN3_HAS_ZLIB
    indexpath += ".gz";
#endif
    std::ofstream ofs{indexpath, std::ios::binary};
    if (ofs)
    {
        // Write the index to disk, including a version string.
#ifdef SEQAN3_HAS_ZLIB
        seqan3::contrib::gz_ostream gzstream(ofs);
        cereal::BinaryOutputArchive oarchive{gzstream};
#else
        cereal::BinaryOutputArchive oarchive{ofs};
#endif
        std::string const version{"1 mars bi_fm_index rna5 collection\n"};
        oarchive(version);
        oarchive(index);
#ifdef SEQAN3_HAS_ZLIB
        gzstream.flush();
#endif
    }
    ofs.close();
}

bool BiDirectionalSearch::read_index(std::filesystem::path const & indexpath)
{
    bool success = false;
#ifdef SEQAN3_HAS_ZLIB
    std::filesystem::path gzindexpath = indexpath;
    gzindexpath += ".gz";
    if (std::filesystem::exists(gzindexpath))
    {
        std::ifstream ifs{gzindexpath, std::ios::binary};
        if (ifs.good())
        {
            seqan3::contrib::gz_istream gzstream(ifs);
            cereal::BinaryInputArchive iarchive{gzstream};
            std::string version;
            iarchive(version);
            assert(version[0] == '1');
            iarchive(index);
            success = true;
        }
        ifs.close();
    }
#endif
    if (!success && std::filesystem::exists(indexpath))
    {
        std::ifstream ifs{indexpath, std::ios::binary};
        if (ifs.good())
        {
            cereal::BinaryInputArchive iarchive{ifs};
            std::string version;
            iarchive(version);
            assert(version[0] == '1');
            iarchive(index);
            success = true;
        }
        ifs.close();
    }
    return success;
}

void BiDirectionalSearch::create_index(std::filesystem::path const & filepath)
{
    if (filepath.empty())
        return;

    std::filesystem::path indexpath = filepath;
    indexpath += ".marsindex";

    // Check whether an index already exists.
    if (read_index(indexpath))
        return;

    // No index found: read genome and create an index.
    if (std::filesystem::exists(filepath))
    {
        // Generate the BiFM index.
        auto seqs = read_genome(filepath);
        index = Index{seqs};
        write_index(indexpath);
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
    assert(!index.empty());
    auto res = seqan3::search(queries.back(), index);
    matches.clear();
    for (auto & result : res)
        matches.emplace_back(result.reference_id(), result.reference_begin_position());
    return matches.size();
}

} // namespace mars
