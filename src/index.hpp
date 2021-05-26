#pragma once

#include <tuple>
#include <vector>

#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>

#include "bi_alphabet.hpp"
#include "index_io.hpp"
#include "motif.hpp"

namespace mars
{

struct Hit
{
    size_t pos;
    uint8_t midx;
    float score;

    Hit(size_t pos, uint8_t midx, float score) : pos{pos}, midx{midx}, score{score} {}
};

//! \brief Provides a bi-directional search step-by-step with backtracking.
class BiDirectionalIndex
{
private:
    //! \brief The index in which the search is performed.
    Index index;

    //! \brief The names of the sequences in the index.
    std::vector<std::string> names;

    //! \brief The history of cursors (needed for backtracking).
    std::vector<seqan3::bi_fm_index_cursor<Index>> cursors;

    //! \brief The history of scores;
    std::vector<float> scores;

    //! \brief The maximum possible motif offset.
    size_t max_offset;

    //! \brief The xdrop parameter.
    unsigned char const xdrop_dist;

public:
    /*!
     * \brief Constructor for a bi-directional search.
     * \param xdrop The xdrop parameter.
     */
    explicit BiDirectionalIndex(unsigned char xdrop):
        index{},
        names{},
        cursors{},
        scores{},
        max_offset{},
        xdrop_dist{xdrop}
    {
        scores.emplace_back(0);
    }

    /*!
     * \brief Create an index of a genome from the specified file.
     * \param filepath The filepath to the file.
     * \return a valid index of the genome.
     * \throws seqan3::file_open_error if neither `filepath` nor `filepath.marsindex` exist.
     *
     * \details
     *
     * This function has two modes:
     *
     * 1. If `filepath.marsindex` exists: Read the already created index from this file.
     * 2. Else if `filepath` exists: Read sequences from this file, create an index
     *    and write the index to `filepath.marsindex`.
     */
    void create(std::filesystem::path const & filepath);

    /*!
     * \brief Append a character to the 5' (left) side of the query.
     * \param item The character to be added.
     * \param left Whether the loop is at the 5' side.
     * \returns whether the operation was successful.
     */
    bool append_loop(std::pair<float, seqan3::rna4> item, bool left);

    /*!
     * \brief Append a character pair at both sides of the query.
     * \param stem_item The score and characters to be added.
     * \returns whether the operation was successful.
     */
    bool append_stem(std::pair<float, bi_alphabet<seqan3::rna4>> stem_item);

    //! \brief Revert the previous append step, which shrinks the query by one or two characters.
    void backtrack();

    /*!
     * \brief Whether the search should be aborted through the xdrop condition.
     * \return True if the score dropped over the previous x elements, false otherwise.
     */
    [[nodiscard]] bool xdrop() const;

    /*!
     * \brief Perform the search with the current query and store the result in `hits`.
     * \param[out] hits The result vector.
     * \param[in] motif The motif for which the results are reported.
     */
    void compute_hits(std::vector<std::vector<Hit>> & hits, StemloopMotif const & motif) const;

    /*!
     * \brief Access a sequence name.
     * \param idx The position of the sequence.
     * \return the name of the queried sequence.
     */
    std::string const & get_name(size_t idx) const
    {
        return names[idx];
    }

    /*!
     * \brief Access the number of sequences in the index.
     * \return the number of sequences
     */
    size_t number_of_seq() const
    {
        return names.size();
    }

    /*!
     * \brief Update the maximum offset of the motif.
     * \param new_offset The proposed new offset value that is used if it is larger.
     */
    void update_max_offset(size_t new_offset)
    {
        if (new_offset > max_offset)
            max_offset = new_offset;
    }

    /*!
     * \brief Retrieve the maximum positional offset of the motif.
     * \return the maximum offset of the motif.
     */
    size_t get_max_offset() const
    {
        return max_offset;
    }
};

} // namespace mars
