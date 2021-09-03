/*
 =============================================================================
 Name        : repeat-coverage-single-cell.cpp
 Author      : Francesca Nadalin
 Version     : 
 Description : Calculate the coverage of repeat elements for each cell barcode
 =============================================================================
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <algorithm>

#include <sys/time.h>
#include <sys/resource.h>

#include "Sam.h"
#include "Bed.h"

#define RESERVE_BLOCK 10000
#define MIN_READ_COUNT 2 // to consider a detected UMI

size_t END_DIFF = 10;

int main(int argc, char * argv[]) {
	
	clock_t time1 = clock();

	if (argc < 5)
	{
		std::cout << "Usage: <bedfile> <bamfile> <outdir> <CB>" << std::endl;
		std::cout << "<bedfile>   contains on the 4th column the repeat element, on the 5th the repeat family, on the 6th the strand" << std::endl;
		std::cout << "<bamfile>   should contain sam fields CB (corrected cellular barcode) and UB (corrected molecular barcode)" << std::endl;
		std::cout << "<outdir>    output directory containing the sparse matrix (matrix, colnames, rownames)" << std::endl;
		std::cout << "<CB>        list of CB to consider, one per line" << std::endl;
		std::cout << "N.B.: both files must be sorted by coordinate and 0-based" << std::endl;
		std::cout << std::endl;
		return 1;
	}

	std::string bed_file(argv[1]);
	std::string bam_file(argv[2]);
	std::string out_dir(argv[3]);
	std::string cb_list_file(argv[4]);

	std::ifstream bed;	
	bed.open(bed_file.c_str());
	if (not bed.is_open())
	{
		std::cerr << "Error while opening file " << bed_file << " for reading" << std::endl;
		return 2;
	}

	std::map<std::string, std::vector<Bed> > beds; // key: reference, value: vector of BED entries with reference = key
	typedef std::map<std::string, std::map<std::string, size_t> > umi_list_t; // key: CB, value = map of UMI and read count
	std::map<std::string, umi_list_t> cb; // key: repeat element, value = umi_list_t
	std::map<std::string, std::string > repeat_families; // key: repeat element, value = repeat family/class the key belongs to

	std::cout << "Reading BED file" << std::endl;

	size_t n = 0;
	for (std::string line; std::getline(bed, line); ) 
	{
		Bed B(line);
		if (beds.count(B.ref()) == 0)
		{
			std::vector<Bed> v; v.reserve(RESERVE_BLOCK);
			beds.insert(std::pair<std::string, std::vector<Bed> >(B.ref(), v));
		}
		if (beds.find(B.ref())->second.capacity() == beds.find(B.ref())->second.size()) 
		{
			beds.find(B.ref())->second.reserve(beds.find(B.ref())->second.capacity() + RESERVE_BLOCK);
		}
		if (cb.count(B.repeat_element()) == 0) 
		{
			umi_list_t l;
			cb.insert(std::pair<std::string, umi_list_t>(B.repeat_element(), l));
			repeat_families.insert(std::pair<std::string, std::string>(B.repeat_element(), B.repeat_family()));
		}
		beds.find(B.ref())->second.push_back(B);
		++n;
		if (n%100000 == 0)
		{
			std::cerr << ".";
			std::cerr.flush();
			if (n%10000000 == 0)
			{
				std::cerr << std::endl;
			}
		}
	}
	bed.close();
	
	std::cout << std::endl << "Done." << std::endl << std::endl;

	std::ifstream cb_list;
	cb_list.open(cb_list_file);
	if (not cb_list.is_open())
	{
		std::cerr << "Error while opening file " << cb_list_file << " for reading" << std::endl;
		return 2;
	}

	std::map<std::string, size_t> CB_map;
	std::vector<std::string> CB;
	CB.reserve(RESERVE_BLOCK);
	n = 0;
	for (std::string c; std::getline(cb_list, c); ) 
	{
		if (CB.capacity() == CB.size()) {
			CB.reserve(CB.capacity() + RESERVE_BLOCK);
		}
		CB.push_back(c);
		CB_map.insert(std::pair<std::string, size_t>(c,n));
		++n;
	}
	CB.resize(CB.size());
	cb_list.close();

	Sam S(bam_file);
	if (not S.is_open())
	{
		std::cerr << "Error while opening file " << bam_file << " for reading" << std::endl;
		return 2;
	}	

	std::cout << "Reading SAM file and collecting info" << std::endl;

	std::map< std::string, std::vector<Bed> >::iterator b = beds.begin();
	size_t i = 0;
	n = 0;
	while (S.read())
	{
		if (S.is_mapped()) 
		{
//			std::cout << "reading sam..." << std::endl;

			// identify a read position with its midpoint computed from the read length
			// => it is independent of the right end of the alignment 
			// => read ordering by coordinate remains consistent, assuming that all the reads have the same length
			if (S.end() - (S.start() + S.length()) < END_DIFF or (S.start() + S.length()) - S.end() < END_DIFF) 
			{
				// center of the read
				size_t mid = S.start() + ((size_t) S.length()/2);
				std::string cb_tag = S.tag_value("CB");
				std::string ub_tag = S.tag_value("UB");
			
//				std::cout << std::endl;
//				std::cout << "ref: " << S.ref() << "	start: " << S.start() << "	length: " << S.length() << "	mid: " << mid << std::endl;
				if (not (cb_tag.empty() or ub_tag.empty()) and CB_map.count(cb_tag) > 0)
				{
//					std::cout << "CB: " << S.tag_value("CB") << "	UB: " << S.tag_value("UB") << std::endl;

					// find the reference of the read in the bam file
					if (S.ref().compare(b->first) != 0)
					{
						if (beds.count(S.ref()) > 0)
						{
							b = beds.find(S.ref());
						}
						i = 0;
					}
//					std::cout << "bed ref: " << b->first << std::endl;
					if (S.ref().compare(b->first) == 0)
					{
						// find the bed entry that contains the read
						while (i < b->second.size() and b->second.at(i).end() <= mid) 
						{
//							std::cout << "bed end: " << b->second.at(i).end() << std::endl;
							++i;
						}
						// the midpoint of the read is inside the repeat and they map to the same strand
						if (i < b->second.size() and b->second.at(i).start() <= mid and
						   not (b->second.at(i).is_reverse() xor S.is_reverse()) and
						   cb.count(b->second.at(i).repeat_element()) > 0)
						{
							// bed interval found!!
//							std::cout << "bed start: " << b->second.at(i).start() << std::endl;
//							std::cout << "bed end: " << b->second.at(i).end() << std::endl;
//							std::cout << "bed interval found!!" << std::endl;
//							std::cout << "repeat element: " << rep_el->first << std::endl;
//							std::cout << "number of CB: " << rep_el->second.size() << std::endl;
//							std::cout << "current CB tag: " << cb_tag << std::endl;
							// find the entry corresponding to the repeat element
							auto rep_el = cb.find(b->second.at(i).repeat_element());
							// if no read from the repeat element is found for this CB, then add a new CB
							if (rep_el->second.count(cb_tag) == 0) 
							{
								std::map<std::string, size_t> m;
								rep_el->second.insert(std::pair<std::string, std::map<std::string, size_t> >(cb_tag, m));
							}
							// find the entry corresponding to the CB
							auto cb_el = rep_el->second.find(cb_tag);
//							std::cout << "CB: " << cb_el->first << std::endl;
							// if no read from the repeat element is found for this (CB,UB) pair, then add a new UB
							if (cb_el->second.count(ub_tag) == 0) 
							{
								cb_el->second.insert(std::pair<std::string, size_t>(S.tag_value("UB"), 0));
							}
							// update the read count
							++(cb_el->second.find(ub_tag)->second);
//							std::cout << "UMI count for (repeat,CB): " << cb_el->second.size() << std::endl;
						}
					}
				}
			}
		}
		++n;
		if (n%100000 == 0)
		{
			std::cerr << ".";
			std::cerr.flush();
			if (n%10000000 == 0)
			{
				std::cerr << std::endl;
			}
		}
	}

	std::cout << std::endl << "Done." << std::endl << std::endl;

	std::cout << "Print to file" << std::endl;

	typedef std::map<std::string, size_t> umi_count_t; // key: CB, value = UMI count
	std::map<std::string, umi_count_t> cb_count; // key: repeat element, value = umi_count_t for each CB
	for (auto ii = cb.begin(); ii != cb.end(); ++ii) 
	{
		umi_count_t umi_count;
		// add a new entry for the repeat element
		cb_count.insert(std::pair<std::string, umi_count_t>(ii->first, umi_count));
		for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj) 
		{
			// count the number of distinct UMI with at least MIN_READ_COUNT reads
			size_t sum = 0;
			for (auto kk = jj->second.begin(); kk != jj->second.end(); ++kk)
			{
				sum += (kk->second >= MIN_READ_COUNT);
			}
			// add a new entry for the CB for the current repeat elements
			// this new entry contains the UMI count
			auto i = cb_count.find(ii->first);
			if (sum > 0)
			{
				i->second.insert(std::pair<std::string, size_t>(jj->first, sum)); 
//				std::cout << "repeat: " << ii->first << "	CB: " << jj->first << "	UMIcount: " << sum << std::endl; 
			}
		}
	}

	// output in sparse matrix format

	// row names = repeat elements
	// col names = CB
	// matrix = UMI count
	std::string out_rows_name(out_dir);
	out_rows_name.append("/features.tsv");
	std::ofstream out_rows;
	out_rows.open(out_rows_name.c_str());
	if (not out_rows.is_open())
	{
		std::cerr << "Error while opening file " << out_rows_name << " for writing" << std::endl;
		return 2;
	}
	std::string out_matrix_name(out_dir);
	out_matrix_name.append("/matrix.mtx");
	std::ofstream out_matrix;
	out_matrix.open(out_matrix_name.c_str());
	if (not out_matrix.is_open())
	{
		std::cerr << "Error while opening file " << out_matrix_name << " for writing" << std::endl;
		return 2;
	}
	// generate the sparse matrix output
	size_t m = 0; // row index
	for (auto i = cb_count.begin(); i != cb_count.end(); ++i)
	{
		// NB: the CB are ordered because the map is a red-black tree (in-order visit)
		out_rows << repeat_families.find(i->first)->second 
			 << "	" << i->first 
			 << "	Gene expression" << std::endl;
		n = 0; // column index
		// scan the CB in the same order as the list in input
		for (auto j = CB.begin(); j != CB.end(); ++j)
		{
			if (i->second.count(*j) > 0)
			{
				size_t c = i->second.find(*j)->second; // entry val
				out_matrix << m << " " << n << " " << c << std::endl;
			}
			++n;
		}
		++m;
	}
	out_rows.close();
	out_matrix.close();
	
	std::cout << "Done." << std::endl << std::endl;

	double time2 = clock();	
	std::cerr << std::endl;
	std::cerr << "Wall time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << std::endl;

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	std::cerr << "RSS: " << usage.ru_maxrss << std::endl;

	return 0;
}

