#include <chrono>
#include <future>
#include <iostream>
#include <vector>

#include "index.hpp"
#include "motif.hpp"
#include "search.hpp"
#include "settings.hpp"

int main(int argc, char ** argv)
{
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    // Set the output stream
    std::ostream out{std::cout.rdbuf()};

    // Parse arguments
    mars::Settings settings{};
    if (!settings.parse_arguments(argc, argv, out))
        return EXIT_FAILURE;

    // Start reading the genome and creating the index asyncronously
    mars::BiDirectionalIndex bds{settings.xdrop};
    std::future<void> index_future = std::async(std::launch::async, &mars::BiDirectionalIndex::create, &bds,
                                                settings.genome_file);

    // Generate motifs from the MSA
    std::vector<mars::StemloopMotif> motifs = mars::create_motifs(settings.alignment_file, settings.threads);

    // Wait for index creation process
    index_future.wait();

    if (!motifs.empty() && !settings.genome_file.empty())
    {
        mars::SearchGenerator search{bds, motifs.front().depth};
        search.find_motifs(motifs);
    }
    else if (mars::verbose > 0)
    {
        std::cerr << "Search step skipped, because there are no motifs, or no genome is provided." << std::endl;
    }

    if (mars::verbose > 0)
    {
        auto const & sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count();
        std::cerr << argv[0] << " has finished after " << sec << " seconds." << std::endl;
    }

    return 0;
}
