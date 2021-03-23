#include <string>                // strings

#include "cli_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct argument_parser_test : public cli_test {};

TEST_F(argument_parser_test, no_options)
{
    cli_test_result result = execute_app("mars");
    std::string expected
    {
        "mars - Motif-based aligned RNA searcher\n"
        "=======================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, "");
}

TEST_F(argument_parser_test, fail_positional_argument)
{
    cli_test_result result = execute_app("mars", data("out.aln"));
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, "Parsing error. Too many arguments provided. Please see -h/--help for more information.\n");
}

//TEST_F(argument_parser_test, with_alignment)
//{
//    cli_test_result result = execute_app("mars", "-a", data("tRNA.aln"));
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
//    EXPECT_EQ(result.err, std::string{});
//}
//
//TEST_F(argument_parser_test, with_argument_verbose)
//{
//    cli_test_result result = execute_app("mars", data("tRNA.aln"), "-v");
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, "> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n");
//    EXPECT_EQ(result.err, "Conversion was a success. Congrats!\n");
//}
//
//TEST_F(argument_parser_test, with_out_file)
//{
//    cli_test_result result = execute_app("mars", data("tRNA.aln"), "-o", "out.fasta");
//    seqan3::sequence_file_input fin{"out.fasta", seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
//
//    // create records to compare
//    using record_type = typename decltype(fin)::record_type;
//    using seqan3::operator""_dna5;
//    std::vector<record_type> records{};
//    records.emplace_back("ACGTTTGATTCGCG"_dna5, std::string{"seq1"});
//    records.emplace_back("TCGGGGGATTCGCG"_dna5, std::string{"seq2"});
//
//    EXPECT_RANGE_EQ(fin, records);
//    EXPECT_EQ(result.exit_code, 0);
//    EXPECT_EQ(result.out, std::string{});
//    EXPECT_EQ(result.err, std::string{});
//}
