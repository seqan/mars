#pragma once

#include <seqan3/std/filesystem>
#include <istream>

#include "multiple_alignment.hpp"

namespace mars
{

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param stream The input stream where the alignment is parsed from.
 * \return The alignment.
 */
msa_type read_msa(std::istream & stream);

/*!
 * \brief Read a CLUSTAL file (*.aln) into a multiple alignment representation.
 * \param filepath The file where the alignment is stored.
 * \return The alignment.
 */
msa_type read_msa(std::filesystem::path const & filepath);

} // namespace mars
