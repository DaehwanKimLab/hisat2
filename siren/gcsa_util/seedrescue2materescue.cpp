/**
 * Seed rescue
 */

#include <utility>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <cstdlib> // exit()
#include "Tools.h"

#define L 10000 // Text row length
//#define DEBUG 1 // Debug output

using namespace std;


void printUsage(char const *name)
{
    cerr << "usage: " << name << " [al_1] [pat_1] [ref]" << endl << endl
         << "Reformats the alignment file to be used for materescue." << endl
         << "\t[al_1]\tOutput file of gcsa_alignment (alignments of single-ends)" << endl
         << "\t[pat_1]\tPattern file (the single-end reads)" << endl
         << "\t[ref]\tThe reference sequence (as a FASTA file)" << endl
         << "Output: Input to materescue.cpp." << endl << endl
         << "Check README for more information." << endl;
}

string readRef(FILE *f)
{
    string ref;
    char row[L];
    while(fgets(row, L, f))
    {
        if (row[0] == '>')
        {
            //cerr << "discarded row: " << row;
            continue;
        }
        if (row[strlen(row)-1] == '\n')
            row[strlen(row)-1] = 0;
        ref.append(row);
    }
    return ref;
}

void revstr(char * t)
{
    char c;
    unsigned n = strlen(t);
    for (unsigned i = 0; i < n / 2; ++i) {
        c = t[i];
        t[i] = t[n - i - 1];
        t[n - i - 1] = c;
    }
}

void complement(char *t)
{
    for (unsigned i = 0; i < strlen(t); ++i)
    {
        switch (t[i])
        {
        case('T'):
            t[i] = 'A';
            break;
        case('G'):
            t[i] = 'C';
            break;
        case('C'):
            t[i] = 'G';
            break;
        case('A'):
            t[i] = 'T';
            break;
        }
    }
}

void readNextRead(FILE *f, int &id, char *pattern, vector<int> &pos, char *row)
{
    if (row[0] == 0)
    {
        id = 1 << 30;
        pattern[0] = 0;
        pos.clear();
        return;
    }
    int p = 0;
    if (sscanf(row, "%d\t%s\t%d\n", &id, pattern, &p) != 3)
    {
        cerr << "error: unable to decode row: " << row;
        exit(1);
    }
    assert(id > 0);
    assert(p > 0);
    pos.clear();
    pos.push_back(p);

    fgets(row, L, f);
    if (feof(f))
    {
        row[0] = 0;
        return;
    }
    int i = -1;
    char tmp[L];
    if (sscanf(row, "%d\t%s\t%d\n", &i, tmp, &p) != 3)
    {
        cerr << "error: unable to decode row: " << row;
        exit(1);
    }
    while (!feof(f) && i == id)
    {
        assert(p > 0);
        pos.push_back(p);        

        fgets(row, L, f);
        if (feof(f))
        {
            row[0] = 0;
            return;
        }
        if (sscanf(row, "%d\t%s\t%u\n", &i, tmp, &p) != 3)
        {
            cerr << "error: unable to decode row: " << row;
            exit(1);
        }
    }
}

void myassert(char *p, char *t)
{
    char tmp[L];
    strcpy(tmp, p);
    tmp[strlen(t)] = 0;
    if (strcmp(tmp,t) == 0)
        return;
    revstr(tmp);
    complement(tmp);
    if (strcmp(tmp,t) != 0)
    {
        cerr << "myassert() failed p and t do not match" << endl;
        cerr << p << endl;
        cerr << t << endl;
        exit(2);
    }
}


bool rescue(string ref, char *pattern, int pn, vector<int> const &pos)
{
    bool res = false;
    MyersEditDistance ed((uchar *)pattern, strlen(pattern), strlen(pattern));
    for(unsigned i = 0; i < pos.size(); ++i)
    {
        // Check both sides of each position
        int p = pos[i];
        int start = 0;
        if (p > (int)strlen(pattern))
            start = p - strlen(pattern);
        string tmp = ref.substr(start, 3*strlen(pattern));
        pair<unsigned,unsigned> r = ed.prefixDist((uchar *)tmp.c_str(), tmp.size());
        if (r.first < 4)
        {
            printf("%d\t%s\t%lu\n", pn, pattern, start + (int)r.second - strlen(pattern));
            res=true;
        }
    }
    return res;
}

int main(int argc, char **argv) 
{
    if (argc != 4)
    {
        printUsage(argv[0]);
        return 1;
    }

    // Init
    FILE *isf = fopen(argv[1], "r");
    if (isf == 0)
    {
        cerr << "error: al_1 file " << argv[1] << " not found" << endl;
        return 1;
    }
    FILE *patternf = fopen(argv[2], "r");
    if (patternf == 0)
    {
        cerr << "error: pat_1 file " << argv[2] << " not found" << endl;
        return 1;
    }
    FILE *reff = fopen(argv[3], "r");
    if (reff == 0)
    {
        cerr << "error: reference file " << argv[3] << " not found" << endl;
        return 1;
    }
    string ref = readRef(reff);
    //cerr << "ref: \"" << ref << "\"" << endl;;
    fclose(reff);

    int f = -1;    // current read id
    char readf[L]; // current read pattern
    vector<int> posf; // vector of current read positions
    char cachef[L]; // first row of the next pattern
    fgets(cachef, L, isf);
    readNextRead(isf, f, readf, posf, cachef);

    char patf[L];
    int nextOutput = 1;
    unsigned rescued = 0, failed = 0;

    while (!feof(isf))
    {        
        // Discard unmapped reads
        while (nextOutput < f)
        {
            fgets(patf, L, patternf); // Skip in pattern file
            ++nextOutput;
        }

        fgets(patf, L, patternf); // Keep pattern files in sync
        patf[strlen(patf)-1]=0; 

#ifdef DEBUG
        cerr << "id = " << f << ", read = " << readf << ", pos =";
        for (vector<int>::iterator it = posf.begin(); it != posf.end(); ++it)
            cerr << ", " << *it;
        cerr << endl;
#endif
        assert(posf.size() > 0);
        myassert(patf, readf);


//        if (posf.size() == 1)
//            cout << f << "\t0\tgi|224589812|ref|nc_000020.10|\t" << posf[0] << "\t255\t*\t*\t0\t0\t*\t*\n"; 
//        else
        {            
            bool res = rescue(ref, patf, f, posf);
            revstr(patf); complement(patf);
            res |= rescue(ref, patf, f, posf);
            revstr(patf); complement(patf);
            if (res)
                ++rescued;
            else
                ++failed;
        }

        readNextRead(isf, f, readf, posf, cachef);
        nextOutput ++;
        if (nextOutput % 100000 == 0)
            cerr << "nextOutput = " << nextOutput << ", " << rescued << " reads rescued with <4e (" << 100*rescued/(rescued+failed) << "%)." << endl;
    }

    cerr << "All input files exausted. " << rescued << " reads rescued." << endl;
    return 0;
}
