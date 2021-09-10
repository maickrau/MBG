#include <vector>
#include <string>
#include <cxxopts.hpp>
#include "MBG.h"

int main(int argc, char** argv)
{
	std::cerr << "MBG " << VERSION << std::endl;
	cxxopts::Options options { "MBG" };
	options.add_options()
		("h,help", "Print help")
		("v,version", "Print version")
		("i,in", "Input reads. Multiple files can be input with -i file1.fa -i file2.fa etc (required)", cxxopts::value<std::vector<std::string>>())
		("o,out", "Output graph (required)", cxxopts::value<std::string>())
		("t", "Number of threads", cxxopts::value<size_t>()->default_value("1"))
		("k", "K-mer size. Must be odd and >=31 (required)", cxxopts::value<size_t>())
		("w", "Window size. Must be 1 <= w <= k-30 (default: k-30)", cxxopts::value<size_t>())
		("a,kmer-abundance", "Minimum k-mer abundance", cxxopts::value<size_t>()->default_value("1"))
		("u,unitig-abundance", "Minimum average unitig abundance", cxxopts::value<double>()->default_value("2"))
		("error-masking", "Error masking", cxxopts::value<std::string>()->default_value("hpc"))
		("include-end-kmers", "Force k-mers at read ends to be included")
		("output-sequence-paths", "Output the paths of the input sequences to a file (.gaf)", cxxopts::value<std::string>())
	;
	auto params = options.parse(argc, argv);
	if (params.count("v") == 1)
	{
		std::cerr << "Version: " << VERSION << std::endl;
		exit(0);
	}
	if (params.count("h") == 1)
	{
		std::cerr << options.help();
		std::cerr << "Options for --error-masking:" << std::endl;
		std::cerr << "\tno\tNo error masking" << std::endl;
		std::cerr << "\thpc\tMask homopolymer errors (default)" << std::endl;
		std::cerr << "\tdinuc\tMask homopolymer and dinucleotide errors" << std::endl;
		std::cerr << "\tmsat\tMask homopolymer and microsatellite errors up to 6bp" << std::endl;
		std::cerr << "\tcollapse\tCollapse homopolymers" << std::endl;
		std::cerr << "\tcollapse-dinuc\tCollapse homopolymers and mask dinucleotide errors" << std::endl;
		std::cerr << "\tcollapse-msat\tCollapse homopolymers and mask microsatellite errors up to 6bp" << std::endl;
		exit(0);
	}
	bool paramError = false;
	if (params.count("i") == 0)
	{
		std::cerr << "Select input reads -i" << std::endl;
		paramError = true;
	}
	if (params.count("o") == 0)
	{
		std::cerr << "Select output filename -o" << std::endl;
		paramError = true;
	}
	if (params.count("k") == 0)
	{
		std::cerr << "Select k-mer size -k" << std::endl;
		paramError = true;
	}
	if (paramError) std::exit(1);
	std::vector<std::string> inputReads = params["i"].as<std::vector<std::string>>();
	std::string outputGraph = params["o"].as<std::string>();
	size_t kmerSize = params["k"].as<size_t>();
	size_t windowSize = 1;
	if (params.count("w") == 1)
	{
		windowSize = params["w"].as<size_t>();
	}
	else
	{
		windowSize = kmerSize - 30;
	}
	size_t minCoverage = params["a"].as<size_t>();
	double minUnitigCoverage = params["u"].as<double>();
	size_t numThreads = params["t"].as<size_t>();
	std::string outputSequencePaths = "";
	ErrorMasking errorMasking = ErrorMasking::Hpc;
	bool includeEndKmers = false;
	std::string errorMaskingStr = "hpc";
	if (params.count("error-masking") == 1)
	{
		errorMaskingStr = params["error-masking"].as<std::string>();
		if (params["error-masking"].as<std::string>() == "no")
		{
			errorMasking = ErrorMasking::No;
		}
		else if (params["error-masking"].as<std::string>() == "hpc")
		{
			errorMasking = ErrorMasking::Hpc;
		}
		else if (params["error-masking"].as<std::string>() == "dinuc")
		{
			errorMasking = ErrorMasking::Dinuc;
		}
		else if (params["error-masking"].as<std::string>() == "collapse")
		{
			errorMasking = ErrorMasking::Collapse;
		}
		else if (params["error-masking"].as<std::string>() == "msat")
		{
			errorMasking = ErrorMasking::Microsatellite;
		}
		else if (params["error-masking"].as<std::string>() == "collapse-dinuc")
		{
			errorMasking = ErrorMasking::CollapseDinuc;
		}
		else if (params["error-masking"].as<std::string>() == "collapse-msat")
		{
			errorMasking = ErrorMasking::CollapseMicrosatellite;
		}
		else
		{
			std::cerr << "unknown parameter for --error-masking: \"" << params["error-masking"].as<std::string>() << "\"" << std::endl;
			paramError = true;
		}
	}
	if (params.count("include-end-kmers") == 1) includeEndKmers = true;
	if (params.count("output-sequence-paths") == 1) outputSequencePaths = params["output-sequence-paths"].as<std::string>();

	if (numThreads == 0)
	{
		std::cerr << "Number of threads cannot be 0" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize > kmerSize)
	{
		std::cerr << "Window size cannot be greater than k-mer size" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize == 0)
	{
		std::cerr << "Window size must be >0" << std::endl;
		paramError = true;
	}
	if (kmerSize % 2 == 0)
	{
		std::cerr << "K-mer size must be odd" << std::endl;
		paramError = true;
	}
	if (kmerSize < 31)
	{
		std::cerr << "Minimum k-mer size is 31" << std::endl;
		paramError = true;
	}
	if (params.count("w") == 1 && windowSize > kmerSize - 30)
	{
		std::cerr << "Window size must be <= k-30. With current k (" << kmerSize << ") maximum w is " << kmerSize - 30 << std::endl;
		paramError = true;
	}
	if (paramError) std::abort();
	
	std::cerr << "Parameters: ";
	std::cerr << "k=" << kmerSize << ",";
	std::cerr << "w=" << windowSize << ",";
	std::cerr << "a=" << minCoverage << ",";
	std::cerr << "u=" << minUnitigCoverage << ",";
	std::cerr << "t=" << numThreads << ",";
	std::cerr << "errormasking=" << errorMaskingStr << ",";
	std::cerr << "endkmers=" << (includeEndKmers ? "yes" : "no");
	std::cerr << std::endl;

	runMBG(inputReads, outputGraph, kmerSize, windowSize, minCoverage, minUnitigCoverage, errorMasking, numThreads, includeEndKmers, outputSequencePaths);
}