// ------------------------------------------------------------------------------------------------------------
// This is MaRs, Motif-based aligned RNA searcher.
// Copyright (c) 2020-2022 Jörg Winkler & Knut Reinert @ Freie Universität Berlin & MPI für molekulare Genetik.
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at https://github.com/seqan/mars.
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <list>
#include <string>
#include <tuple>
#include <vector>

// The submodule lib/ipknot has no namespace

/*!
 * \brief Compute the secondary structure of a given multiple structural alignment (MSA).
 * \param names The IDs of the MSA.
 * \param seqs The sequences of the MSA.
 * \return two vectors which hold the base pairs and pseudoknot levels.
 */
std::pair<std::vector<int>, std::vector<int>> run_ipknot(std::list<std::string> const & names,
                                                         std::list<std::string> const & seqs);
