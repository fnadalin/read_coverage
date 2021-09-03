#ifndef SAM_H
#define SAM_H

#define JUNCTION_COORD_TOL 5 // bp around the splicing junction for a spanning read
#define MIN_INTRON_SPAN 3 // min bp in the intron for an unspliced read to be considered non-spanning

#include "samtools/sam.h"
#include "common.h"

typedef std::pair <size_t, size_t> junction; // donor start (0-based) ; acceptor end (0-based, past-the-end)

#include <iostream>
#include <string>
class Sam {

	public:

	const size_t FLAG_UNMAP = 4;

	Sam ();
	Sam (const std::string & s);
	~Sam();

	void open (const std::string & s);
	bool is_open () { return _is_open; }
	bool is_sam () { return _is_sam; }
	bool read ();
	std::string id ();
	std::string sequence (bool strand);
	bool is_first ();
	bool is_reverse ();
	// bool is_strand_ok (bool pair_end); // reverse complement if the sequence is not compatible with a paired read
	std::string ref () const;
	size_t start () const; // 0-based
	size_t length () const;
	size_t end () const; // 0-based, past-the-end
	bool is_mapped () const;
	std::vector<junction> junctions () const;
	std::string tag_value (const char tag[2]) const;

	private:

	samfile_t * _filesam;
	bam1_t * _b;
	unsigned _is_open: 1, _is_sam: 1;

};

#endif 
