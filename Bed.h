#ifndef BED_H
#define BED_H

#include "common.h"

#include <string>
#include <sstream>


class Bed {

	public:

	Bed (const std::string & s);
	~Bed();

	void read (const std::string & s);
	std::string ref () const { return _ref; }
	size_t start () const { return _start; }
	size_t end () const { return _end; }
	std::string name () const { return _name; }
	std::string score () const { return _score; }
	std::string repeat_element () const { return name(); } // alias for name()
	std::string repeat_family () const { return score(); } // alias for score()
	bool is_reverse () const { return _strand == '-'; }

	private:

	std::string   _ref;
	size_t        _start;
	size_t        _end;
	std::string   _name;
	std::string   _score;
	char          _strand;
	size_t        _thickStart;
	size_t        _thickEnd;
	std::string   _itemRgb;
	size_t        _blockCount;
	size_t        _blockSizes;
	size_t        _blockStarts;
	
	size_t        _num_blocks;

};

#endif 
