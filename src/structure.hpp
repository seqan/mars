#pragma once

#include <list>
#include <string>
#include <tuple>
#include <vector>

#include "multiple_alignment.hpp"

// The submodule lib/ipknot has no namespace

/*!
 * \brief Compute the secondary structure of a given multiple structural alignment (MSA).
 * \param names The IDs of the MSA.
 * \param seqs The sequences of the MSA.
 * \return two vectors which hold the base pairs and pseudoknot levels.
 */
std::pair<std::vector<int>, std::vector<int>> run_ipknot(std::list<std::string> const & names,
                                                         std::list<std::string> const & seqs);

namespace mars
{

/*!
 * \brief Compute the secondary structure of a given multiple structural alignment.
 * \param msa The multiple structural alignment.
 */
void compute_structure(Msa & msa);

}
