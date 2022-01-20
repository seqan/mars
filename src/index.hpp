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
     * \brief Query the number of sequences that have been stored in the index.
     * \return the number of sequences.
     */
    size_t seq_count() const
    {
        return names.size();
    }

    /*!
     * \brief Access a sequence name.
     * \param idx The position of the sequence.
     * \return the name of the `idx`th sequence.
     */
    std::string const & seq_name(size_t idx) const
    {
        return names[idx];
    }

    Index const & raw() const
    {
        return index;
    }
};

} // namespace mars
