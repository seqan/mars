// IPknot interface

#pragma once

#include <list>
#include <string>
#include <tuple>
#include <vector>

std::pair<std::vector<int>, std::vector<int>> run_ipknot(std::list<std::string> const & names,
                                                         std::list<std::string> const & seqs);
