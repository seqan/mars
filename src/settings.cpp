#include <thread>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "settings.hpp"

namespace mars
{

unsigned short verbose{0};

bool Settings::parse_arguments(int argc, char ** argv, std::ostream & out)
{
    seqan3::argument_parser parser{"mars", argc, argv};
    parser.info.short_description = "Motif-based aligned RNA searcher";
    parser.info.author = "Jörg Winkler";
    parser.info.version = "1.0.0";
    parser.info.date = "March 2021";
    parser.info.examples.emplace_back("./mars structural_rna.aln -g genome.fasta -o out.txt");
    parser.info.description.emplace_back("MaRs is a tool that reads a structural multiple RNA alignment "
                                         "(e.g. from LaRA) and derives fuzzy stem loop descriptors from it. "
                                         "These descriptors are then subject to a search in an indexed database or "
                                         "sequence and MaRs returns the hits where the RNA structure is found, "
                                         "accompanied with a quality value for each hit.");

    // Parser
    parser.add_option(alignment_file, 'a', "alignment",
                      "Alignment file of structurally aligned RNA sequences.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"msa", "aln", "fasta", "fa", "sth", "stk"}});

    parser.add_option(genome_file, 'g', "genome",
                      "A sequence file containing one or more sequences.");

    //output path as option, otherwise output is printed
    parser.add_option(result_file, 'o', "output",
                      "The output file for the results. If empty we print to stdout.");

    parser.add_option(min_score_per_motif, 's', "scorefilter",
                      "Minimum score per motif that a hit must achieve. Influences the output of low-scoring hits.");

    parser.add_option(xdrop, 'x', "xdrop",
                      "The xdrop parameter. Smaller values increase speed but we will find less matches.");

    parser.add_option(threads, 'j', "threads",
                      "Use the number of specified threads. Value 0 tries to detect the maximum number.");

    parser.add_option(verbose, 'v', "verbose",
                      "Level of printing status information.");

    try
    {
        parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                    // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return false;
    }

    if (threads == 0u)
    {
        unsigned int nthreads = std::thread::hardware_concurrency();
        threads = nthreads != 0u ? nthreads : 1u;
    }
    if (verbose > 0)
        std::cerr << "Number of threads: " << threads << std::endl;

    file_stream.rdbuf()->open(result_file, std::ios_base::out);
    if (!result_file.empty())
        out.rdbuf(file_stream.rdbuf());

    if (!out)
    {
        seqan3::debug_stream << "Failed to open the output file " << result_file << "\n";
        return false;
    }

    return true;
}

} // namespace mars
