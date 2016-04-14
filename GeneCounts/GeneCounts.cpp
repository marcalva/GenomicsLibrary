//
//  GeneCounts.cpp
//
//  Created by Marcus Alvarez on 01/12/16.
//  Pajukanta Lab
//

#include "GeneCounts.h"

std::vector<string> parse_string(const std::string& in, const char& delim, std::vector<std::string>& out)
{
	out.clear();
    std::stringstream ss(in);
    std::string token;

    while (std::getline(ss, token, delim))
	    out.push_back(token);

    return out;
}

Gene::Gene(std::string Name)
{
	Name = Name;
	Count = 0;
}

bool GeneCounts::ReadBED(std::string BEDFileName)
{
	_bedFileName = BEDFileName;
	Tabix bedFs(_bedFileName);
	
	vector<string> bed_vector;
	string bed_line;
	while ( bedFs.getNextLine(bed_line) )
	{
		parse_string(bed_line, '\t', bed_vector);
		string name = bed_vector.at(3);
		Gene g(name);
		_countsTable.insert(pair<std::string,Gene>(name,g));
	}
	return true;
}

bool GeneCounts::ReadGTF(std::string GTFFileName)
{
	return true;
}

bool GeneCounts::OpenBAM(std::string BAMFileName)
{
	_bamfs.Open(BAMFileName);
	
	if ( !_bamfs.IsOpen() )
	{
		cout << "Could not open BAM file " << BAMFileName << endl;
		return false;
	}
	
	string bami = BAMFileName + ".bai";
	if ( ! _bamfs.OpenIndex(bami) )
	{
		cout << "Creating index for " << BAMFileName << endl;
		_bamfs.CreateIndex();
	}
	
	if ( ! _bamfs.HasIndex() )
	{
		cout << "No access to BAM index" << endl;
		return false;
	}
	return true;
}

void GeneCounts::CountReads()
{	
	Tabix bedFs(_bedFileName);
	
	BamTools::BamAlignment al;
	while (_bamfs.GetNextAlignment(al))
	{
		string featureName;
		bool featureDefined = false;
		
		string alignChr = _bamfs.GetReferenceData().at(al.RefID).RefName;
		int alignStart = al.Position;
		string readName = al.Name;
		
		bool countRead = true;
	
		// For each CIGAR alignment
		for (unsigned int i = 0; i < al.CigarData.size(); i++)
		{
			char t = al.CigarData.at(i).Type;
			unsigned int l = al.CigarData.at(i).Length;
		
			// hard or soft clip, or insertion is not part of alignment so pass
			if (t == 'H' || t == 'S' || t == 'I')
				continue;
			else if ( t == 'D' || t == 'N' )
				alignStart += l;
			else if (t == 'M')
			{
				int alignEnd = alignStart + l; // Half open-ended
				alignStart++;
				alignEnd++;
				Interval readInterval(alignChr, alignStart, alignEnd);
				string start = to_string(alignStart);
				string end = to_string(alignEnd);
				string region = alignChr + ":" + start + "-" + end;
				
				// Check Overlap
				GenomicInterval allIntervals;
				list<Interval> overlapIntervals;
				vector<string> bed_vector;
				string bed_line;
				bedFs.setRegion(region);
				while ( bedFs.getNextLine(bed_line) )
				{
					parse_string(bed_line, '\t', bed_vector);
					Interval ti(bed_vector.at(0), atoi(bed_vector.at(1).c_str()), atoi(bed_vector.at(2).c_str()), bed_vector.at(3));
					allIntervals.addInterval(ti);
				}
				overlapIntervals = allIntervals.intersect(readInterval);
				if ( overlapIntervals.size() == 0 ) { countRead = false; }
				else if ( overlapIntervals.size() == 1 )
				{
					if ( !featureDefined ) { featureName = overlapIntervals.front().Name(); }
					else if ( featureDefined && featureName != overlapIntervals.front().Name() ) { countRead = false; }
				}
				else if ( overlapIntervals.size() > 1 ) // Check if all overlapping intervals have the same feature name
				{
					list<Interval>::const_iterator oii = overlapIntervals.begin();
					string n = oii->Name();
					oii++;
					bool skip = false;
					for (; oii != overlapIntervals.end(); ++oii)
					{
						if ( oii->Name() != n ) { skip=true;}
					}
					if (skip) {countRead = false;}
					else
					{
						if ( !featureDefined ) { featureName = n; }
						else if (n != featureName) { countRead = false; }
					}
				}
				alignStart += l;
			} // end if CIGAR == 'M'
			else
			{
				cout << "Unrecognized CIGAR string: " << t << endl;
				exit(1);
			}
		}
		
		if (countRead)
		{
			set<string>::const_iterator it;
			it = _countsTable[featureName].ReadNames.find(readName);
			if (it == _countsTable[featureName].ReadNames.end())
			{
				_countsTable[featureName].ReadNames.insert(readName);
				_countsTable[featureName].Count++;
			}
			else
			{
				_countsTable[featureName].ReadNames.erase(it);
				_countsTable[featureName].Count++;
			}
		}
	}
}

bool GeneCounts::WriteTable(std::string OutFileName)
{
	return true;
}