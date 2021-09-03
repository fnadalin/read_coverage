#include "Bed.h"

Bed::Bed (const std::string & s)
{
	std::vector<std::string> v;
	v.reserve(12); // max number of fields in BED
	std::istringstream is(s);
	std::string el;
	char delim = '	';
	while (std::getline(is, el, delim))
	{
		v.push_back(el);
	}

	_ref = "0";
	_start = 0;
	_end = 0;
	_name = ".";
	_score = ".";

	size_t n = 0;
	if (v.size() > n) _ref = v.at(n++);
	if (v.size() > n) _start = atoi(v.at(n++).c_str());
	if (v.size() > n) _end = atoi(v.at(n++).c_str());
	if (v.size() > n) _name = v.at(n++);
	if (v.size() > n) _score = v.at(n++);
	if (v.size() > n) _strand = v.at(n++)[0];
	if (v.size() > n) _thickStart = atoi(v.at(n++).c_str());
	if (v.size() > n) _thickEnd = atoi(v.at(n++).c_str());
//	if (v.size() > n) _itemRgb = v.at(n++);
//	if (v.size() > n) _blockCount = atoi(v.at(n++).c_str());
//	if (v.size() > n) _blockSizes = atoi(v.at(n++).c_str());
//	if (v.size() > n) _blockStarts = atoi(v.at(n++).c_str());

}

// do nothing
Bed::~Bed()
{

}


