/*
 ============================================================================
 Name        : read-coverage.cpp
 Author      : Francesca Nadalin
 Version     : 
 Description : Calculate the distribution of reads on genomic locations
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

	if (argc < 4)
	{
		std::cout << "Usage: <bedfile> <bamfile> <outfile>" << std::endl;
		std::cout << "N.B.: both files must be sorted by coordinate and 0-based" << std::endl;
		std::cout << "Output file: BED format with number of reads in the 4th column" << std::endl;
		return 1;
	}

	std::string bed_file(argv[1]);
	std::string bam_file(argv[2]);
	std::string out_file(argv[3]);

	std::ifstream bed;	
	bed.open(bed_file.c_str());
	if (not bed.is_open())
	{
		std::cerr << "Error while opening file " << bed_file << " for reading" << std::endl;
	}

	std::map<std::string, std::vector<Bed> > beds; // key: reference, value: vector of BED entries with reference = key
	std::map<std::string, std::vector<size_t> > read_count; // key: reference, value: vector of the number of reads whose center falls in the BED entry pointed by key in beds
	for (std::string line; std::getline(bed, line); ) 
	{
		Bed B(line);
		if (beds.count(B.ref()) == 0)
		{
			std::vector<Bed> v; v.reserve(RESERVE_BLOCK);
			std::vector<size_t> w; w.reserve(RESERVE_BLOCK);
			beds.insert(std::pair<std::string, std::vector<Bed> >(B.ref(), v));
			read_count.insert(std::pair<std::string, std::vector<size_t> >(B.ref(), w));
		}
		if (beds.find(B.ref())->second.capacity() == beds.find(B.ref())->second.size()) 
		{
			beds.find(B.ref())->second.reserve(beds.find(B.ref())->second.capacity() + RESERVE_BLOCK);
			read_count.find(B.ref())->second.reserve(read_count.find(B.ref())->second.capacity() + RESERVE_BLOCK);
		}
		beds.find(B.ref())->second.push_back(B);
		read_count.find(B.ref())->second.push_back(0);
	}
	bed.close();

	Sam S(bam_file);
	if (not S.is_open())
	{
		std::cerr << "Error while opening file " << bam_file << " for reading" << std::endl;
	}	

	std::map< std::string, std::vector<Bed> >::iterator b = beds.begin();
	std::map< std::string, std::vector<size_t> >::iterator r = read_count.begin();
	size_t i = 0;
	size_t n = 0;
	while (S.read())
	{
		if (S.is_mapped()) 
		{
//			std::cout << "reading sam..." << std::endl;
			// center of the read
			size_t mid = S.start() + ((size_t) S.length()/2);
		
//			std::cout << "ref: " << S.ref() << "	start: " << S.start() << "	length: " << S.length() << "	mid: " << mid << std::endl;
			// find the reference of the read in the bam file
			if (S.ref().compare(b->first) != 0)
			{
				if (beds.count(S.ref()) > 0)
				{
					b = beds.find(S.ref());
					r = read_count.find(S.ref());
				}
				i = 0;
			}
//			std::cout << "bed ref: " << b->first << std::endl;
			if (S.ref().compare(b->first) == 0)
			{
				// find the bed entry that contains the read
				while (i < b->second.size() and b->second.at(i).end() <= mid) 
				{
//					std::cout << "bed end: " << b->second.at(i).end() << std::endl;
					++i;
				}
				if (i < b->second.size() and b->second.at(i).start() <= mid)
				{
					// bed interval found!!
//					std::cout << "bed end: " << b->second.at(i).end() << std::endl;
//					std::cout << "bed start: " << b->second.at(i).start() << std::endl;
//					std::cout << "bed interval found!!" << std::endl;
					++(r->second.at(i));
				
				}
			}
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

	// print
	std::ofstream out;
	out.open(out_file.c_str());
	if (not out.is_open())
	{
		std::cerr << "Error while opening file " << out_file << " for writing" << std::endl;
	}
	b = beds.begin();
	r = read_count.begin();
	while (b != beds.end() and r != read_count.end())
	{
		for (size_t i = 0; i < b->second.size(); ++i)
		{
			out << b->second.at(i).ref() << "	" 
				<< b->second.at(i).start() << "	" 
				<< b->second.at(i).end() << "	" 
				<< b->second.at(i).name() << "	" 
				<< b->second.at(i).score() << "	"
				<< r->second.at(i) << std::endl;
		}
		++b;
		++r;
	}
	out.close();
	

	double time2 = clock();	
	std::cout << std::endl;
	std::cout << "Wall time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << std::endl;

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	std::cout << "RSS: " << usage.ru_maxrss << std::endl;

	return 0;
}

