#include <seqan3/argument_parser/all.hpp>

#include "settings.hpp"

namespace mars
{

std::unique_ptr<thread_pool::ThreadPool> pool;
std::mutex mutex_cerr;
Settings settings{};

bool Settings::parse_arguments(int argc, char ** argv)
{
    unsigned int nthreads{std::thread::hardware_concurrency()};
    nthreads = nthreads != 0u ? nthreads : 1u;

    seqan3::argument_parser parser{"mars", argc, argv};
    parser.info.version = "1.0.0";
    parser.info.short_description = "Motif-based aligned RNA searcher";
    parser.info.author = "Jörg Winkler";
    parser.info.email = "j.winkler@fu-berlin.de";
    parser.info.date = "April 2022";
    parser.info.url = "https://github.com/seqan/mars";
    parser.info.short_copyright = "2020-2022 Jörg Winkler and Knut Reinert, FU-Berlin; "
                                  "released under the 3-clause BSD license.";
    parser.info.long_copyright = "BSD 3-Clause License\n"
                                 "\n"
                                 "Copyright (c) 2020-2022, Jörg Winkler and Knut Reinert,\n"
                                 "Freie Universität Berlin & MPI für molekulare Genetik.\n"
                                 "All rights reserved.\n"
                                 "\n"
                                 "Redistribution and use in source and binary forms, with or without\n"
                                 "modification, are permitted provided that the following conditions are met:\n"
                                 "\n"
                                 "1. Redistributions of source code must retain the above copyright notice, this\n"
                                 "   list of conditions and the following disclaimer.\n"
                                 "\n"
                                 "2. Redistributions in binary form must reproduce the above copyright notice,\n"
                                 "   this list of conditions and the following disclaimer in the documentation\n"
                                 "   and/or other materials provided with the distribution.\n"
                                 "\n"
                                 "3. Neither the name of the copyright holder nor the names of its\n"
                                 "   contributors may be used to endorse or promote products derived from\n"
                                 "   this software without specific prior written permission.\n"
                                 "\n"
                                 "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n"
                                 "AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n"
                                 "IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\n"
                                 "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE\n"
                                 "FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\n"
                                 "DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\n"
                                 "SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
                                 "CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,\n"
                                 "OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
                                 "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";
    parser.info.description.emplace_back("MaRs is a tool that reads a structural multiple RNA alignment "
                                         "(e.g. from LaRA) and derives fuzzy stem loop descriptors from it. "
                                         "These descriptors are then subject to a search in an indexed database or "
                                         "sequence and MaRs returns the hits where the RNA structure is found, "
                                         "accompanied with a quality value for each hit.");
    parser.info.synopsis.emplace_back("./mars structuralRNA.aln -g genome.fasta -o out.txt");

    // Parser
    parser.add_subsection("Input data:");
    parser.add_option(genome_file, 'g', "genome",
                      "A sequence file containing one or more sequences.");

#if SEQAN3_WITH_CEREAL
    parser.add_option(alignment_file, 'a', "alignment",
                      "Alignment file of structurally aligned RNA sequences, "
                      "or a motif file to restore previously calculated motifs.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"msa", "aln", "sth", "stk", "sto", "mmo"}});
#else
    parser.add_option(alignment_file, 'a', "alignment",
                      "Alignment file of structurally aligned RNA sequences.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"msa", "aln", "sth", "stk", "sto"}});
#endif

    //output path as option, otherwise output is printed
    parser.add_subsection("Output options:");
    parser.add_option(result_file, 'o', "output",
                      "The output file for the results. If empty we print to stdout.");

#if SEQAN3_WITH_CEREAL
    parser.add_option(motif_file, 'm', "motif", "File for storing the motifs.", seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"mmo"}});
#endif

    parser.add_option(structator_file, 'r', "rssp", "Output rssp file for the Structator program.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"pat"}});

    parser.add_option(min_score_per_motif, 's', "scorefilter",
                      "Minimum score per motif that a hit must achieve. Influences the output of low-scoring hits.");

    parser.add_option(evalue_filter, 'e', "evalue",
                      "Filter based on e-value instead of score.");

    parser.add_option(verbose, 'v', "verbose",
                      "Level of printing status information.");

    parser.add_subsection("Performance options:");
    parser.add_option(prune, 'p', "prune",
                      "Prune the search if occurence is lower than p% of expected.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0,100});

    parser.add_option(xdrop, 'x', "xdrop",
                      "The xdrop parameter. Smaller values increase speed but we will find less matches.");

    parser.add_option(exterior, 'l', "exterior",
                      "Consider long exterior and multibranch loops as well.");

#ifdef SEQAN3_HAS_ZLIB
    parser.add_option(compress_index, 'z', "gzip",
                      "Use gzip compression for the index file.");
#endif

    parser.add_option(nthreads, 'j', "threads",
                      "Use the number of specified threads.");

    try
    {
        parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                    // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return false;
    }

    pool = std::make_unique<thread_pool::ThreadPool>(nthreads);
    return true;
}

} // namespace mars
