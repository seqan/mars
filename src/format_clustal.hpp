#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/utility/views/type_reduce.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "multiple_alignment.hpp"

namespace mars
{

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the sequences.
 * \param stream The input stream where the alignment is parsed from.
 * \return The alignment.
 */
template<seqan3::alphabet alphabet_type>
MultipleAlignment<alphabet_type> read_clustal_file(std::istream & stream)
{
    MultipleAlignment<alphabet_type> msa;

    // Define a lambda function to decide if a character is legal.
    auto check_legal_alphabet = [] (char const c)
    {
        using legal_alphabet_type = std::conditional_t<seqan3::nucleotide_alphabet<alphabet_type>,
                                                       seqan3::rna15,
                                                       seqan3::aa27>;
        auto constexpr is_legal_alph = seqan3::char_is_valid_for<seqan3::gapped<legal_alphabet_type>>;
        if (!is_legal_alph(c))
            throw seqan3::parse_error{"Encountered an unexpected letter: char_is_valid_for<" +
                                      seqan3::detail::type_name_as_string<legal_alphabet_type> +
                                      "> evaluated to false on " + seqan3::detail::make_printable(c)};
        return c;
    };

    auto stream_view = seqan3::detail::istreambuf(stream);

    // skip initial whitespace and check if file starts with "CLUSTAL"
    seqan3::detail::consume(stream_view | seqan3::detail::take_until(!seqan3::is_space));
    if (!std::ranges::equal(stream_view | seqan3::detail::take_exactly_or_throw(7), std::string{"CLUSTAL"}))
        throw seqan3::parse_error{"Expected to read 'CLUSTAL' in the beginning of the file."};

    // skip rest of the line and go to the first block
    seqan3::detail::consume(stream_view | seqan3::detail::take_line);
    seqan3::detail::consume(stream_view | seqan3::detail::take_until(!seqan3::is_space));

    std::size_t idx = 0;
    bool first_block = true;

    do // read block-wise
    {
        // parse the sequence name
        std::string name{};
        std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_blank),
                          std::cpp20::back_inserter(name));

        if (!msa.names.empty() && name == msa.names.front())
        {
            // if the first name is found again, initiate the next block
            first_block = false;
            idx = 0;
        }

        if (first_block)
        {
            // add new name and sequence
            msa.names.push_back(name);
            msa.sequences.push_back(std::vector<seqan3::gapped<alphabet_type>>{});
        }
        else
        {
            // check for inconsistencies
            assert(msa.names.size() > idx);
            if (!std::ranges::equal(name, msa.names[idx])) // validate the sequence name
                throw seqan3::parse_error{"Expected to read '" + msa.names[idx] + "' in the input file."};
        }
        // go to the beginning of sequence
        seqan3::detail::consume(stream_view | seqan3::detail::take_until_or_throw(!seqan3::is_blank));

        // copy the sequence
        std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_space)
                          | std::views::filter(!seqan3::is_digit)
                          | std::views::transform(check_legal_alphabet)
                          | seqan3::views::char_to<seqan3::gapped<alphabet_type>>,
                          std::cpp20::back_inserter(msa.sequences[idx]));

        // go to the next line
        seqan3::detail::consume(stream_view | seqan3::detail::take_line);

        // consume line if it starts with whitespace
        if (seqan3::is_space(*stream_view.begin()))
            seqan3::detail::consume(stream_view | seqan3::detail::take_line);

        // move to the next sequence name or find EOF
        seqan3::detail::consume(stream_view | seqan3::detail::take_until(seqan3::is_eof || (!seqan3::is_space)));
        ++idx;
    } while (std::istreambuf_iterator<char>{stream} != std::istreambuf_iterator<char>{});

    return std::move(msa);
}

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the sequences.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
template<seqan3::alphabet alphabet_type>
MultipleAlignment<alphabet_type> read_clustal_file(std::filesystem::path const & filepath)
{
    // Open filepath as stream.
    std::ifstream stream(filepath, std::ios_base::in | std::ios::binary);
    if (!stream.good())
        throw seqan3::file_open_error{"Could not open file " + filepath.string() + " for reading."};

    auto result = read_clustal_file<alphabet_type>(stream);
    stream.close();
    return std::move(result);
}

} // namespace mars
