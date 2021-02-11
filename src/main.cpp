#include <chrono>
#include <vector>

#include "index.hpp"
#include "input_output.hpp"
#include "motif.hpp"
#include "search.hpp"
#include "settings.hpp"
#include "structure.hpp"

int main(int argc, char ** argv)
{
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // Set the output stream
    std::ostream out{std::cout.rdbuf()};

    // Parse arguments
    mars::Settings settings{};
    if (!settings.parse_arguments(argc, argv, out))
        return EXIT_FAILURE;

    // Read the alignment
    mars::Msa msa = mars::read_msa(settings.alignment_file);

    // Compute an alignment structure
    auto structure = mars::compute_structure(msa);

    // Find the stem loops
    std::vector<mars::StemloopMotif> motifs = mars::detect_stemloops(structure.first, structure.second);

    // Create a structure motif for each stemloop
    std::for_each(motifs.begin(), motifs.end(), [&msa, &structure] (mars::StemloopMotif & motif)
    {
        motif.analyze(msa, structure.first);
    });

    // Read genome or quit
    if (settings.genome_file.empty())
    {
        for (auto & motif : motifs)
            out << motif;
    }
    else
    {
        mars::Index index = mars::create_index(settings.genome_file);
        mars::SeqNum const depth = msa.sequences.size();
        mars::SearchGenerator search{std::move(index), depth, settings.xdrop};
        for (auto & motif : motifs)
        {
            std::cerr << motif;
            search.find_motif(motif);
        }
    }

    auto const & sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count();
    std::cerr << "Run time " << sec << " seconds.\n";

    return 0;
}
