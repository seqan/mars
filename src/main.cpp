#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include "index.hpp"
#include "motif.hpp"
#include "search.hpp"
#include "settings.hpp"

int main(int argc, char ** argv)
{
    std::chrono::steady_clock::time_point tm0 = std::chrono::steady_clock::now();

    // Parse arguments
    if (!mars::settings.parse_arguments(argc, argv))
        return EXIT_FAILURE;

    // Start reading the genome and creating the index asyncronously
    mars::BiDirectionalIndex bds{};
    auto future_index = mars::pool->submit(&mars::BiDirectionalIndex::create, &bds);

    // Generate motifs from the MSA
    std::vector<mars::StemloopMotif> motifs = mars::create_motifs();
    auto future_json = mars::pool->submit(mars::store_motifs, motifs);
    auto future_rssp = mars::pool->submit(mars::store_rssp, motifs);

    // Wait for index creation process
    future_index.wait();

    if (!motifs.empty() && !mars::settings.genome_file.empty())
    {
        mars::SearchGenerator search{bds};
        search.find_motifs(motifs);

        auto print_results = [&bds, &search] (std::ostream & out)
        {
            if (!search.get_locations().empty())
                out << " " << std::left << std::setw(35) << "sequence name" << "\t" << "index" << "\t"
                    << "pos" << "\t" << "n" << "\t" << "score" << std::endl;
            for (mars::MotifLocation const & loc : search.get_locations())
                out << ">" << std::left << std::setw(35) << bds.get_name(loc.sequence) << "\t" << loc.sequence << "\t"
                    << loc.position << "\t" << +loc.num_stemloops << "\t" << loc.score << std::endl;
        };

        if (!mars::settings.result_file.empty())
        {
            logger(1, "Writing results ==> " << mars::settings.result_file << std::endl);
            std::ofstream file_stream(mars::settings.result_file);
            print_results(file_stream);
            file_stream.close();
        }
        else
        {
            print_results(std::cout);
        }
    }
    else if (motifs.empty())
    {
        logger(1, "There are no motifs: skipping search step." << std::endl);
    }
    else
    {
        logger(1, "No genome provided: skipping search step." << std::endl);
    }
    future_json.wait();
    future_rssp.wait();

    // print run time
    auto const sec = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - tm0).count();
    logger(1, argv[0] << " has finished after " << sec << " seconds." << std::endl);

    return 0;
}
