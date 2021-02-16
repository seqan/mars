#include <chrono>
#include <future>
#include <vector>

#ifdef MARS_WITH_OPENMP
    #include <omp.h>
#endif

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

    std::future<mars::Index> index_future = std::async(std::launch::async, mars::create_index, settings.genome_file);

    // Read the alignment
    mars::Msa msa = mars::read_msa(settings.alignment_file);

    // Compute an alignment structure
    auto structure = mars::compute_structure(msa);

    // Find the stem loops
    std::vector<mars::StemloopMotif> motifs = mars::detect_stemloops(structure.first, structure.second);

    // Create a structure motif for each stemloop
    #pragma omp parallel for num_threads(2)
    for (size_t idx = 0; idx < motifs.size(); ++idx)
        motifs[idx].analyze(msa, structure.first);

    // Read genome or quit
    mars::Index index = index_future.get();
    if (index.empty())
    {
        for (auto & motif : motifs)
            out << motif;
    }
    else
    {
        mars::SeqNum const depth = msa.sequences.size();
        mars::SearchGenerator search{std::move(index), depth, settings.xdrop};
        search.find_motifs(motifs);
    }

    auto const & sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count();
    std::cerr << "Run time " << sec << " seconds.\n";

    return 0;
}
