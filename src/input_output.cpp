#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "format_clustal.hpp"
#include "input_output.hpp"

namespace mars
{

msa_type read_msa(std::istream & stream)
{
    return read_clustal_file<seqan3::rna15>(stream);
}

msa_type read_msa(std::filesystem::path const & filepath)
{
    return read_clustal_file<seqan3::rna15>(filepath);
}

index_type read_genome(std::filesystem::path const & filepath)
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
        // Read rna5 sequences.
        struct my_traits : seqan3::sequence_file_input_default_traits_dna
        {
            using sequence_alphabet = seqan3::rna5;
            using sequence_legal_alphabet = seqan3::rna15;
        };

        std::vector<seqan3::rna5_vector> seqs{};
        for (auto && record : seqan3::sequence_file_input<my_traits>{filepath})
            seqs.push_back(std::move(seqan3::get<seqan3::field::seq>(record)));

        // Generate the BiFM index.
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

} // namespace mars
