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
    parser.add_option(genome_file, 'g', "genome",
                      "A sequence file containing one or more sequences.");

#if SEQAN3_WITH_CEREAL
    parser.add_option(alignment_file, 'a', "alignment",
                      "Alignment file of structurally aligned RNA sequences, "
                      "or a motif file to restore previously calculated motifs.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"msa", "aln", "sth", "stk", "sto", "json"}});

    parser.add_option(motif_file, 'm', "motif", "File for storing the motifs.", seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"json"}});
#else
    parser.add_option(alignment_file, 'a', "alignment",
                      "Alignment file of structurally aligned RNA sequences.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"msa", "aln", "sth", "stk", "sto"}});
#endif

    //output path as option, otherwise output is printed
    parser.add_option(result_file, 'o', "output",
                      "The output file for the results. If empty we print to stdout.");

    parser.add_option(structator_file, 'r', "rssp", "Output file for rssp output for the Structator program.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"pat"}});

    parser.add_option(prune, 'p', "prune",
                      "Prune the search if occurence is lower than p% of expected.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{0,100});

    parser.add_option(max_evalue, 'e', "evalue",
                      "Maximum e-value for result list. Influences the output of low-scoring hits.");

    parser.add_option(xdrop, 'x', "xdrop",
                      "The xdrop parameter. Smaller values increase speed but we will find less matches.");

    parser.add_option(nthreads, 'j', "threads",
                      "Use the number of specified threads.");

    parser.add_option(verbose, 'v', "verbose",
                      "Level of printing status information.");

#ifdef SEQAN3_HAS_ZLIB
    parser.add_option(compress_index, 'z', "gzip",
                      "Use gzip compression for the index file.");
#endif

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
