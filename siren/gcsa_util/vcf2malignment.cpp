// Convert SNP data to simplified multiple alignment
#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
using namespace std;

#define R 6  // Number of strands in the resulting multiple alignment
#define MAXLENGTH 1000000 // Max length of a line

struct snp {
    unsigned i; // Start and end position in reference (0-based)
    unsigned j;
    string alt; // Alternative strand
    bool outputted;
    int strain;
};

bool is_dna(string const &t)
{
    for (string::const_iterator it = t.begin(); it != t.end(); ++it)
        if (*it != 'A' && *it != 'C' 
            && *it != 'G' && *it != 'T' && *it != '-' && *it != 'N')
            return false;
    return true;
}

void skip_to(unsigned pos, unsigned offset, string const &sequence)
{
    for (; offset < pos; ++offset)
    {
        for (int i = 0; i < R; ++i)
            cout << sequence[offset];
                
        cout << "\n";
    }
}

void flush(unsigned offset, vector<snp> const &output)
{
    unsigned i = 1;
    bool exausted = false;
    while (!exausted)
    {
        exausted = true;
        for (vector<snp>::const_iterator it = output.begin(); it != output.end(); ++it)
        {
            if (it->j-1 == offset && offset - it->i + i < it->alt.size() && it->outputted)
                exausted = false;
        }

        if (!exausted)
        {
            cout << '-';
            int r = 0;
            for (vector<snp>::const_iterator it = output.begin(); it != output.end(); ++it)
            {
                if (it->j-1 == offset && it->strain != -1)
                {
                    assert(it->outputted);
                    assert (it->strain >= r);
                    while (r < it->strain)
                    {
                        ++r;
                        cout << '-';
                    }

                    if (offset - it->i + i < it->alt.size())
                        cout << it->alt[ offset - it->i + i ];
                    else
                        cout << '-';
                    ++r;
                }
            }
            while ((++r) < R)
                cout << '-';
            cout << endl;
        }
        ++i;
    }
}

bool mysort (snp const &i, snp const &j) 
{ 
    return (i.strain < j.strain);
}

unsigned fill(unsigned offset, vector<snp>  &output, string const &sequence)
{
    if (offset < output[0].i)
        skip_to(output[0].i, offset, sequence);

    offset = output[0].i;
    unsigned end = offset;
    for (vector<snp>::const_iterator it = output.begin(); it != output.end(); ++it)
        if (end < it->j)
            end = it->j;

    set<int> freeStrains;
    for (int i = 0; i < R-1; ++i)
        freeStrains.insert(i);

    while (offset < end)
    {
        // Assign strain numbers
        for (vector<snp>::iterator it = output.begin(); it != output.end(); ++it)
        {
            if (it->i == offset) 
            {
                assert (!it->outputted);
                if (freeStrains.empty())
                {
                    cerr << "Hit the R limit at " << offset << endl;
                    continue;
                }
                
                set<int>::iterator free = freeStrains.begin();
                it->strain = *free;
                it->outputted = true;
                freeStrains.erase(free);
            }
        }

        // Output
        cout << sequence[offset];
        int r = 0;
        sort(output.begin(), output.end(), mysort);

        for (vector<snp>::iterator it = output.begin(); it != output.end(); ++it)
        {
            if (it->strain != -1)
            {
                /**/
                assert (it->outputted);
                assert (it->strain >= r);
                while (r < it->strain)
                {
                    ++r;
                    cout << sequence[offset];
                }

                if (offset == it->j)
                    cout << sequence[offset];
                else if (offset - it->i < it->alt.size())
                    cout << it->alt[offset - it->i];
                else
                    cout << '-';
                ++r;
            }
        }

        while ((++r) < R)
            cout << sequence[offset];
        cout << "\n";

        flush(offset, output);

        // Free strains
        for (vector<snp>::iterator it = output.begin(); it != output.end(); ++it)
            if (it->strain != -1 && offset == it->j)
            {
                assert (freeStrains.find(it->strain) == freeStrains.end());
                freeStrains.insert(it->strain); // Free the strain
                it->strain = -1;
            }
        

        ++offset;
    }

    for (vector<snp>::iterator it = output.begin(); it != output.end(); ++it)
        it->strain = -1;
    return offset;
}



int main(int argc, char* argv[]) 
{
    if (argc != 3)
    {
        cerr << "usage: " << argv[0] << " snpfile chromno < fasta > outputfile" << endl;
        return 1;
    }

    string chromno(argv[2]);

    // Read FASTA
    char row[MAXLENGTH];
    fgets(row, MAXLENGTH, stdin);
    if (row[0] != '>')
    {
        cerr << "error: input not in FASTA format?" << endl;
        return 2;
    }
    cerr << "Reading seq: " << row;

    string sequence;
    unsigned size = 0;

    fgets(row, MAXLENGTH, stdin);
    while (!feof(stdin))
    {
        assert(row[0] != '>');
        sequence.append(row, strlen(row)-1);
        size += strlen(row)-1;
        fgets(row, MAXLENGTH, stdin);
    }
    cerr << "Total length = " << sequence.size() << endl;
    assert (size == sequence.size());

    // Read SNPs
    unsigned offset = 0;

    FILE *fp = fopen(argv[1], "rt");
    if (!fp)
    {
        cerr << "error: unable to read input file " << argv[1] << endl;
        return 3;
    }
    fgets(row, MAXLENGTH, fp);
    cerr << "Discarding the first row: " << row << endl;
    fgets(row, MAXLENGTH, fp);

    vector<snp> snps;
    unsigned prevpos = 0;
    while (!feof(fp))
    {
        char chrom[250];
        unsigned pos = 0;
        if (sscanf(row, "%s\t%u", chrom, &pos) != 2)
        {
            cerr << "error: unable to read input file row; " << row << endl;
            return 4;
        }
        if (chromno.compare(chrom) != 0)
        {
            fgets(row, MAXLENGTH, fp);
            continue; // Discard
        }
        assert(pos > 0);
        assert(pos >= prevpos);
        prevpos = pos;
        --pos;
        
        char *tmp = row;
        for (int i = 0; i < 3 && *tmp != 0; ++i, ++tmp) // Skip 3 tabs
        {
            while (*tmp != '\t' && *tmp != 0)
                ++tmp;
            assert (*tmp != 0);
        }

        unsigned len = 0;
        while (tmp[len] != '\t' && tmp[len] != 0) 
            ++len;
        if (tmp[len] == 0)
            cerr << "Invalid row? " << row;
        assert(tmp[len] != 0);
        assert(len >= 1);

        string ref(tmp, len);

        //cerr << "got ref: \"" << ref << "\"" << endl;

        tmp += len+1;
        len = 0;
        while (tmp[len] != '\t' && tmp[len] != 0 && tmp[len] != '\n' && tmp[len] != '\r') 
            ++len;
//        assert(tmp[len] == '\n' || tmp[len] == '\r');

        string alt(tmp, len);
        //cerr << "got alt: \"" << alt << "\"" << endl;

        // Sanity check
        for (unsigned i = 0; i < ref.size(); ++i)
            assert(sequence[pos + i] == ref[i]);

        if (alt.compare("<DEL>") == 0)
            alt = "-";

        if (is_dna(alt))
        {
            snp ns = { pos, pos+ref.size(), alt, false, -1 };
//            cout << "i , j = " << pos << ", " << pos+ref.size() << endl;
            snps.push_back(ns);
        }
        else
        {
            cerr << "Discarding row (is_dna failed): " << row;
        }
        fgets(row, MAXLENGTH, fp);
    }

    cerr << "number of snps = " << snps.size() << endl;
    
    unsigned i = 0;
    offset = 0;
    while (i < snps.size())
    {
        unsigned j = i+1;
        unsigned end = snps[i].j;
        while (snps[j].i < end && j < snps.size())
        {
            if (end < snps[j].j)
                end = snps[j].j;
            ++j;
        }

        // Output from offset to MAX(snps[i'].j)
        vector<snp> output;
        for (unsigned k = 0; k < j-i; ++k)
            output.push_back(snps[i+k]);
        unsigned check = fill(offset, output, sequence);
        for (unsigned k = 0; k < j-i; ++k)
            snps[i+k].outputted = output[k].outputted;

        assert(check == end);
        offset = end; // MAX (snps[i'].j);
        i = j;
    }
    skip_to(sequence.size(), offset, sequence);
    
    cerr << "number of snps = " << snps.size() << endl;
    cerr << "offset = " << offset << endl;
    unsigned total_output = 0, total_filtered = 0;

    for(unsigned i = 0; i < snps.size(); ++i)
        if (snps[i].outputted)
            total_output ++;
        else
            total_filtered ++;
        
    cerr << "outputted  = " << total_output << endl
         << "filtered = " << total_filtered << endl;
    fclose(fp);
    return 0;
}
