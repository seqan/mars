#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>

namespace mars
{

/*!
 * \brief A multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the contained sequences.
 */
template<seqan3::alphabet alphabet_type>
struct multiple_alignment
{
    std::vector<std::vector<seqan3::gapped<alphabet_type>>> sequences;
    std::vector<std::string> names;
};

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \tparam alphabet_type The alphabet type of the sequences.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
template<seqan3::alphabet alphabet_type>
multiple_alignment<alphabet_type> read_clustal_file(std::filesystem::path const & filepath)
{
    multiple_alignment<alphabet_type> msa;

    // Define a lambda function to decide if a character is legal.
    auto check_legal_alphabet = [] (char const c)
    {
        using legal_alphabet_type = std::conditional_t<seqan3::nucleotide_alphabet<alphabet_type>,
                                                       seqan3::rna15,
                                                       seqan3::aa27>;
        auto constexpr is_legal_alph = seqan3::is_in_alphabet<seqan3::gapped<legal_alphabet_type>>;
        if (!is_legal_alph(c))
            throw seqan3::parse_error{"Encountered an unexpected letter: " + is_legal_alph.msg +
                                      " evaluated to false on " + seqan3::detail::make_printable(c)};
        return c;
    };

    // Open filepath as stream.
    std::ifstream stream(filepath, std::ios_base::in | std::ios::binary);
    if (!stream.good())
        throw seqan3::file_open_error{"Could not open file " + filepath.string() + " for reading."};

    auto stream_view = seqan3::views::istreambuf(stream);

    // skip initial whitespace and check if file starts with "CLUSTAL"
    seqan3::detail::consume(stream_view | seqan3::views::take_until(!seqan3::is_space));
    if (!std::ranges::equal(stream_view | seqan3::views::take_exactly_or_throw(7), std::string{"CLUSTAL"}))
        throw seqan3::parse_error{"Expected to read 'CLUSTAL' in the beginning of file " + filepath.generic_string()};

    // skip rest of the line and go to the first block
    seqan3::detail::consume(stream_view | seqan3::views::take_line);
    seqan3::detail::consume(stream_view | seqan3::views::take_until(!seqan3::is_space));

    do // read the first block
    {
        msa.names.emplace_back("");
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_blank),
                          std::cpp20::back_inserter(msa.names.back()));
        seqan3::detail::consume(stream_view | seqan3::views::take_until_or_throw(!seqan3::is_blank));

        msa.sequences.push_back(std::vector<seqan3::gapped<alphabet_type>>{});
        std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_space)
                          | std::views::filter(!seqan3::is_digit)
                          | std::views::transform(check_legal_alphabet)
                          | seqan3::views::char_to<seqan3::gapped<alphabet_type>>,
                          std::cpp20::back_inserter(msa.sequences.back()));

        seqan3::detail::consume(stream_view | seqan3::views::take_line);
    } while (!seqan3::is_space(*stream_view.begin()));

    // move to next block
    seqan3::detail::consume(stream_view | seqan3::views::take_line);
    seqan3::detail::consume(stream_view | seqan3::views::take_until(!seqan3::is_space));

    // read until EOF
    while (std::istreambuf_iterator<char>{stream} != std::istreambuf_iterator<char>{})
    {
        // for each sequence in a block
        for (auto && [name, seq] : seqan3::views::zip(msa.names, msa.sequences))
        {
            // validate the sequence name
            if (!std::ranges::equal(stream_view | seqan3::views::take_exactly_or_throw(name.size()), name))
                throw seqan3::parse_error{"Expected to read '" + name + "' in file " + filepath.generic_string()};

            // skip to sequence and append to alignment
            seqan3::detail::consume(stream_view | seqan3::views::take_until_or_throw(!seqan3::is_blank));
            std::ranges::copy(stream_view | seqan3::views::take_until_or_throw(seqan3::is_space)
                              | std::views::filter(!seqan3::is_digit)
                              | std::views::transform(check_legal_alphabet)
                              | seqan3::views::char_to<seqan3::gapped<alphabet_type>>,
                              std::cpp20::back_inserter(seq));
            seqan3::detail::consume(stream_view | seqan3::views::take_line);
        }

        // move to next block
        seqan3::detail::consume(stream_view | seqan3::views::take_line);
        seqan3::detail::consume(stream_view | seqan3::views::take_until(!seqan3::is_space));
    }

    stream.close();

    return std::move(msa);
}

} // namespace mars
