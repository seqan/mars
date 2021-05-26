#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/gz_istream.hpp>
    #include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

#include "index_io.hpp"

namespace mars
{

void read_genome(std::vector<seqan3::dna4_vector> & seqs,
                 std::vector<std::string> & names,
                 std::filesystem::path const & filepath)
{
    struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
        using sequence_legal_alphabet = seqan3::dna15;
    };
    typedef seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> SeqInput;

    for (auto & [seq, name] : SeqInput{filepath})
    {
        seqs.push_back(std::move(seq));
        names.push_back(std::move(name));
    }
}

void write_index(Index const & index, std::vector<std::string> const & names, std::filesystem::path & indexpath)
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
        std::string const version{"1 mars bi_fm_index<dna4,collection>\n"};
        oarchive(version);
        oarchive(index);
        oarchive(names);
#ifdef SEQAN3_HAS_ZLIB
        gzstream.flush();
#endif
    }
    ofs.close();
}

bool read_index(Index & index, std::vector<std::string> & names, std::filesystem::path & indexpath)
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
            iarchive(names);
            success = true;
            indexpath = gzindexpath;
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
            iarchive(names);
            success = true;
        }
        ifs.close();
    }
    return success;
}

} // namespace mars
