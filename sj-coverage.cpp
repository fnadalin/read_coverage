/*
 ============================================================================
 Name        : read-coverage.cpp
 Author      : Francesca Nadalin
 Version     : 
 Description : Extract the reads spanning splicing junctions
 ============================================================================
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <map>

#include <sys/time.h>
#include <sys/resource.h>

#include "Sam.h"
#include "Bed.h"

#define RESERVE_BLOCK 100

int main(int argc, char * argv[]) {
	
	clock_t time1 = clock();

	if (argc < 5)
	{
		std::cout << "Usage: <junctions.bed> <bamfile> <outfile_spanning> <outfile_nonspanning>" << std::endl;
		std::cout << "<junctions.bed>   BED file with junctions, and with junction ID in the last column" << std::endl;
		std::cout << "Output file: TSV file with spanning read ID and list of junction IDs, comma-separated" << std::endl;
		std::cout << "Output file: TSV file with non-spanning read ID and list of junction IDs, comma-separated" << std::endl;
		return 1;
	}

	std::string bed_file(argv[1]);
	std::string bam_file(argv[2]);
	std::string out_sj_file(argv[3]);
	std::string out_nsj_file(argv[4]);

	std::ifstream bed;	
	bed.open(bed_file.c_str());
	if (not bed.is_open())
	{
		std::cerr << "Error while opening file " << bed_file << " for reading" << std::endl;
		return 2;
	}

	std::vector<Bed> beds; // vector of BED entries
	beds.reserve(RESERVE_BLOCK);
	for (std::string line; std::getline(bed, line); ) 
	{
		Bed B(line);
		if (beds.capacity() == beds.size()) 
		{
			beds.reserve(beds.capacity() + RESERVE_BLOCK);
		}
		beds.push_back(B);
	}
	bed.close();

	Sam S(bam_file);
	if (not S.is_open())
	{
		std::cerr << "Error while opening file " << bam_file << " for reading" << std::endl;
	}	

	std::ofstream out_sj;
	out_sj.open(out_sj_file.c_str());
	if (not out_sj.is_open())
	{
		std::cerr << "Error while opening file " << out_sj_file << " for reading" << std::endl;	
		return 2;
	}

	std::ofstream out_nsj;
	out_nsj.open(out_nsj_file.c_str());
	if (not out_nsj.is_open())
	{
		std::cerr << "Error while opening file " << out_nsj_file << " for reading" << std::endl;	
		return 2;
	}

	int n = 0;
	while (S.read())
	{
		if (S.is_mapped()) 
		{
//			std::cout << "reading sam..." << std::endl;
			// center of the read
//			std::cout << "start: " << S.start() << std::endl;
//			std::cout << "end: " << S.end() << std::endl;
			std::vector<junction> jj = S.junctions();
			size_t block_start = S.start();
			std::string j_list_nspan;
//			std::cout << "junctions:" << std::endl;
			if (jj.size() > 0)
			{
				std::string j_list_span;
				for (auto j : jj) {
//					std::cout << j.first << " " << j.second << std::endl;
					size_t block_end = j.first;
					for (auto b : beds)
					{
						if (S.ref().compare(b.ref()) == 0)
						{
							if (j.first > b.start() - JUNCTION_COORD_TOL and
								j.first < b.start() + JUNCTION_COORD_TOL and
								j.second > b.end() - JUNCTION_COORD_TOL and
								j.second < b.end() + JUNCTION_COORD_TOL)
							{
								// the read overlaps the junction
								// here it is spliced: it encompasses the two spliced exons
								if (j_list_span.size() > 0) j_list_span.append(",");
								j_list_span.append(b.name());
							}
							else 
							{
								if ((block_start < b.start() - JUNCTION_COORD_TOL and
									 block_end > b.start() + JUNCTION_COORD_TOL + MIN_INTRON_SPAN) or 
									(block_start < b.end() - JUNCTION_COORD_TOL - MIN_INTRON_SPAN and
									 block_end > b.start() + JUNCTION_COORD_TOL))
								{
									// the read overlaps the junction
									// here it is not spliced: it spans one exon and the next intron
									if (j_list_nspan.size() > 0) j_list_nspan.append(",");
									j_list_nspan.append(b.name());
								}
							}
						}
					}
					block_start = j.second;
				}
				if (j_list_span.size() > 0) out_sj << j_list_span << "	" << S.id() << std::endl;
			}	
			size_t block_end = S.end();
			for (auto b : beds)
			{
				if (S.ref().compare(b.ref()) == 0)
				{
					if ((block_start < b.start() - JUNCTION_COORD_TOL and
						 block_end > b.start() + JUNCTION_COORD_TOL + MIN_INTRON_SPAN) or 
						(block_start < b.end() - JUNCTION_COORD_TOL - MIN_INTRON_SPAN and
						 block_end > b.start() + JUNCTION_COORD_TOL))
					{
						// the read overlaps the junction
						// here it is not spliced: it spans one exon and the next intron
						if (j_list_nspan.size() > 0) j_list_nspan.append(",");
						j_list_nspan.append(b.name());
					}
				}
			}
			if (j_list_nspan.size() > 0) out_nsj << j_list_nspan << "	" << S.id() << std::endl;
		}
		++n;
		if (n%100000 == 0)
		{
			std::cout << ".";
			std::cout.flush();
			if (n%10000000 == 0)
			{
				std::cout << std::endl;
			}
		}
	}

	return 0;

	double time2 = clock();	
	std::cout << std::endl;
	std::cout << "Wall time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << std::endl;

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	std::cout << "RSS: " << usage.ru_maxrss << std::endl;

	return 0;
}

