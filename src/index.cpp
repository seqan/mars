// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/gz_istream.hpp>
    #include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

#include "index.hpp"
#include "settings.hpp"

namespace mars
{

void BiDirectionalIndex::read_genome(std::vector<seqan3::dna4_vector> & seqs)
{
    struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
    };
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using SeqFileInput = seqan3::sequence_file_input<dna4_traits, fields>;

    auto parse = [this, &seqs] (auto && reader)
    {
        reader.options.truncate_ids = true;

        for (auto & [seq, name] : reader)
        {
            seqs.push_back(std::move(seq));
            names.push_back(std::move(name));
        }
    };

    try
    {
        SeqFileInput reader{settings.genome_file};
        parse(reader);
        return;
    }
    catch (std::invalid_argument & e)
    {
        logger(1, "Could not interpret the file suffix " << settings.genome_file.extension()
                  << ", trying to parse fasta." << std::endl);
    }
    std::ifstream ifs{settings.genome_file};
    if (ifs)
    {
        SeqFileInput reader{ifs, seqan3::format_fasta()};
        parse(reader);
    }
    ifs.close();
}

void BiDirectionalIndex::write_index(std::filesystem::path & indexpath)
{
#ifdef SEQAN3_HAS_ZLIB
    if (settings.compress_index)
    {
        indexpath += ".gz";
        std::ofstream ofs{indexpath, std::ios::binary};
        if (ofs)
        {
            // Write the index to disk, including a version string.
            seqan3::contrib::gz_ostream gzstream(ofs);
            cereal::BinaryOutputArchive oarchive{gzstream};
            std::string const version{"1 mars bi_fm_index<dna4,collection>\n"};
            oarchive(version);
            oarchive(index);
            oarchive(names);
            gzstream.flush();
        }
        ofs.close();
    }
    else
#endif
    {
        std::ofstream ofs{indexpath, std::ios::binary};
        if (ofs)
        {
            // Write the index to disk, including a version string.
            cereal::BinaryOutputArchive oarchive{ofs};
            std::string const version{"1 mars bi_fm_index<dna4,collection>\n"};
            oarchive(version);
            oarchive(index);
            oarchive(names);
        }
        ofs.close();
    }
}

bool BiDirectionalIndex::read_index(std::filesystem::path & indexpath)
{
    bool success = false;
    if (std::filesystem::exists(indexpath))
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
#ifdef SEQAN3_HAS_ZLIB
    if (!success)
    {
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
    }
#endif
    return success;
}

void BiDirectionalIndex::create()
{
    if (settings.genome_file.empty())
        return;

    std::filesystem::path indexpath = settings.genome_file;
    indexpath += ".marsindex";

    // Check whether an index already exists.
    if (read_index(indexpath))
    {
        logger(1, "Using existing index <== " << indexpath << std::endl);
        return;
    }

    // No index found: read genome and create an index.
    if (std::filesystem::exists(settings.genome_file))
    {
        std::vector<seqan3::dna4_vector> seqs{};
        read_genome(seqs);
        logger(1, "Read " << seqs.size() << " genome sequences <== " << settings.genome_file << std::endl);
        if (!seqs.empty())
        {
            // Generate the BiFM index.
            index = Index{seqs};
            write_index(indexpath);
            logger(1, "Created index ==> " << indexpath << std::endl);
        }
    }
    else
    {
        std::ostringstream err_msg{};
        err_msg << "Could not find the genome file <== " << settings.genome_file << "[.marsindex";
#ifdef SEQAN3_HAS_ZLIB
        err_msg << "[.gz]";
#endif
        err_msg << "]";
        throw seqan3::file_open_error(err_msg.str());
    }
}

} // namespace mars
