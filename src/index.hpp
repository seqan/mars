// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/std/filesystem>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

namespace mars
{

//! \brief The type of a bi-directional index over the 4-letter DNA alphabet.
using Index = seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>;

//! \brief Provides a bi-directional search step-by-step with backtracking.
class BiDirectionalIndex
{
private:
    //! \brief The index in which the search is performed.
    Index index;

    //! \brief The names of the sequences in the index.
    std::vector<std::string> names;

    /*!
     * \brief Read a FASTA file of sequences.
     * \param[out] seqs The object where the sequences can be stored.
     */
    void read_genome(std::vector<seqan3::dna4_vector> & seqs);

    /*!
     * \brief Archive an index and store it in a file on disk.
     * \param[in,out] indexpath The path of the index output file.
     */
    void write_index(std::filesystem::path & indexpath);

    /*!
     * \brief Unarchive an index and read it from a file on disk.
     * \param[in] indexpath The path of the index input file.
     * \return whether an index could be parsed.
     */
    bool read_index(std::filesystem::path & indexpath);

public:
    /*!
     * \brief Create an index of a genome.
     * \throws seqan3::file_open_error if neither `genome_file` nor `genome_file.marsindex` exist.
     *
     * \details
     *
     * This function has two modes:
     *
     * 1. If `genome_file.marsindex` exists: Read the already created index from this file.
     * 2. Else if `genome_file` exists: Read sequences from this file, create an index
     *    and write the index to `genome_file.marsindex`.
     */
    void create();

    /*!
     * \brief Access the sequence names.
     * \return the name vector of the sequences.
     */
    std::vector<std::string> const & get_names() const
    {
        return names;
    }

    /*!
     * \brief Access the underlying bi-FM index.
     * \return the raw index (without metadata).
     */
    Index const & raw() const
    {
        return index;
    }
};

} // namespace mars
