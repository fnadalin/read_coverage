#include "Sam.h"

Sam::Sam () {
	_is_open = false;
	_b = new bam1_t;
	_b->data = (uint8_t*)malloc(1024);
}

Sam::Sam (const std::string & s) {

	unsigned short dot = s.find_last_of(".");
	std::string ext = s.substr(dot+1);
	_is_sam = (ext.compare("sam") == 0);
	
	void * aux;
	if (_is_sam){
		std::cout << "is_sam" << std::endl;
	    _filesam = samopen(s.c_str(), "r", aux);
	} else if (ext == "bam") {
	    _filesam = samopen(s.c_str(), "rb", aux);
	}

	_is_open = (_filesam != NULL);
	if (_is_open) {
		_b = new bam1_t;
		_b->data = (uint8_t*)malloc(1024);
	}

	_is_open = (_filesam != NULL);

	std::cout << "l_text: " << _filesam->header->l_text << std::endl;
	std::cout << "n_targets: " << _filesam->header->n_targets << std::endl;
}


Sam::~Sam () {
	if (_is_open) samclose(_filesam);
	if (_b != NULL) {
		free(_b->data);
		delete _b;
	}
}


void Sam::open (const std::string & s) {

	unsigned short dot = s.find_last_of(".");
	if (dot >= s.size()) {
		return;
	}
	std::string ext = s.substr(dot+1);
	_is_sam = (ext.compare("sam") == 0);

	void * aux;
	if (_is_sam){
		_filesam = samopen(s.c_str(), "r", aux);
	} else if (ext == "bam") {
		_filesam = samopen(s.c_str(), "rb", aux);
	}

	_is_open = (_filesam != NULL);
}


bool Sam::read () {
	return (samread (_filesam, _b) > 0);
}


std::string 
Sam::id () {
	std::string id (bam1_qname(_b));
	return id;
}


std::string 
Sam::sequence (bool strand) {

	if (_b == NULL) {
		return std::string();
	}
	uint8_t * seq =  bam1_seq(_b);
	int32_t l = _b->core.l_qseq;
	std::string ts;
	ts.reserve(l+1);
	for (int32_t i = 0; i < l; i++) {
		switch (bam1_seqi(seq, i)) {
		case 1:
			ts.push_back('A');
			break;
		case 2:
			ts.push_back('C');
			break;
		case 4:
			ts.push_back('G');
			break;
		case 8:
			ts.push_back('T');
			break;
		default:
			ts.push_back('N');
			break;
		}
	}
	if (strand) {
		return ts;
	} else {
		return reverse_complement_standalone_str (ts);
	}
}


bool Sam::is_first () {
	return (_b == NULL) ? false : _b->core.flag & BAM_FREAD1;
}


bool Sam::is_reverse () {
		return (_b == NULL) ? false : _b->core.flag & BAM_FREVERSE;
}


/*
bool Sam::is_strand_ok (bool pair_end) {
	bool innie = (is_first () and not is_reverse ()) or (not is_first() and is_reverse());
	return pair_end ? innie : not innie;
}
*/


std::string 
Sam::ref () const
{
	return (_filesam->header->n_targets > _b->core.tid) ? std::string(_filesam->header->target_name[_b->core.tid]) : std::string();
}


size_t
Sam::start () const
{
	return (_b->core.pos);
}


size_t
Sam::length () const
{
	return (_b->core.l_qseq);
}


size_t
Sam::end () const
{
	uint32_t* cigar = bam1_cigar(_b);
	size_t pos = _b->core.pos;
//	std::cout << "cigar: " << cigar << std::endl;
	for (size_t k = 0; k < _b->core.n_cigar; ++k)
	{
		int cop = cigar[k] & BAM_CIGAR_MASK; // operation
		int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
		switch (cop) {
		case BAM_CMATCH:
//			printf("M");
			// printf("[%d-%d]", pos, pos + cl - 1);
			pos += cl;
			break;

		case BAM_CHARD_CLIP:
//			printf("H");
      		/* printf("[%d]", pos);  // No coverage */
      		/* pos is not advanced by this operation */
      		break;

		case BAM_CSOFT_CLIP:
//			printf("S");
			/* printf("[%d]", pos);  // No coverage */
			/* pos is not advanced by this operation */
			break;

		case BAM_CDEL:
//			printf("D");
			/* printf("[%d-%d]", pos, pos + cl - 1);  // Spans positions, No Coverage */
			pos += cl;
			break;

		case BAM_CPAD:
//			printf("P");
			/* printf("[%d-%d]", pos, pos + cl - 1); // Spans positions, No Coverage */
			pos += cl;
			break;

		case BAM_CINS:
//			printf("I");
			/* printf("[%d]", pos); // Special case - adds <cl> bp "throughput", but not genomic position "coverage" */
			/* How you handle this is application dependent */
			/* pos is not advanced by this operation */
			break;

		case BAM_CREF_SKIP:
//			printf("N");
			/* printf("[%d-%d]", pos, pos + cl - 1); /* Spans positions, No Coverage */
			pos += cl;
			break;

		default:
			fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
			printf("?");
		}
	}

	return pos;
}

bool 
Sam::is_mapped () const
{
	return (not (_b->core.flag & FLAG_UNMAP));
}


std::vector <junction> 
Sam::junctions () const 
{
	std::vector <junction> v;
	v.reserve(5);
	uint32_t* cigar = bam1_cigar(_b);
	size_t pos = _b->core.pos;
	for (size_t k = 0; k < _b->core.n_cigar; ++k)
	{
		int cop = cigar[k] & BAM_CIGAR_MASK; // operation
		int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
		switch (cop) {
		case BAM_CMATCH:
			pos += cl;
			break;

		case BAM_CHARD_CLIP:
      		break;

		case BAM_CSOFT_CLIP:
			break;

		case BAM_CDEL:
			pos += cl;
			break;

		case BAM_CPAD:
			pos += cl;
			break;

		case BAM_CINS:
			break;

		case BAM_CREF_SKIP: // intron
			v.push_back(std::make_pair(pos, pos+cl));
			pos += cl;
			break;

		default:
			fprintf(stderr, "Unhandled cigar_op %d:%d\n", cop, cl);
			printf("?");
		}
	}

	return v;
}


std::string
Sam::tag_value (const char tag[2]) const
{
	uint8_t* tag_val_char = bam_aux_get(_b, tag);
	std::string tag_val;
	if (tag_val_char != 0)
		tag_val = (char*) ++tag_val_char; // skip the first character since it contains the middle field
	
	return tag_val;
}

