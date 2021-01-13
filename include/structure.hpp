#pragma once

#include <list>
#include <string>
#include <tuple>
#include <vector>

#include "multiple_alignment.hpp"

// The submodule lib/ipknot has no namespace
std::pair<std::vector<int>, std::vector<int>> run_ipknot(std::list<std::string> const & names,
                                                         std::list<std::string> const & seqs);

namespace mars
{

std::pair<std::vector<int>, std::vector<int>> compute_structure(msa_type const & msa);

}
