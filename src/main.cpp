#include <seqan3/std/filesystem>

#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "input_output.hpp"
#include "motif.hpp"
#include "structure.hpp"

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"mars", argc, argv};
    parser.info.short_description = "Motif-based aligned RNA searcher";
    parser.info.author = "Jörg Winkler";
    parser.info.version = "1.0.0";
    parser.info.date = "September 2020";
    parser.info.examples.emplace_back("./mars structural_rna.aln -g genome.fasta -o out.txt");
    parser.info.description.emplace_back("MaRs is a tool that reads a structural multiple RNA alignment "
                                         "(e.g. from LaRA) and derives fuzzy stem loop descriptors from it. "
                                         "These descriptors are then subject to a search in an indexed database or "
                                         "sequence and MaRs returns the hits where the RNA structure is found, "
                                         "accompanied with a quality value for each hit.");

    // Declarations for argument parser
    std::filesystem::path alignment_file{};
    std::filesystem::path genome_file{};
    std::filesystem::path result_file{};

    // Parser
    parser.add_positional_option(alignment_file, "Alignment file of structurally aligned RNA sequences.",
                                 seqan3::input_file_validator{{"msa", "aln", "fasta", "fa", "sth", "stk"}});
    parser.add_option(genome_file, 'g', "genome", "A sequence file containing one or more sequences.");
    //output path as option, otherwise output is printed
    parser.add_option(result_file, 'o', "output", "The file for the result output, if empty we print to stdout.");

    try
    {
        parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                    // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return -1;
    }

    // Read the alignment
    mars::msa_type msa = mars::read_msa(alignment_file);

    // Compute an alignment structure
    auto && [bpseq, plevel] = mars::compute_structure(msa);

    // Find the stem loops
    mars::stemloop_type stemloops = mars::detect_stem_loops(bpseq, plevel);
    seqan3::debug_stream << stemloops << "\n";

    for (std::pair<int, int> const & pos : stemloops)
        mars::analyze_stem_loop(msa, bpseq, pos);

    return 0;
}
