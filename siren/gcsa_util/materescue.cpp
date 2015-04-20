/**
 * Mate rescue
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
    cerr << "usage: " << name << " [al_1] [al_2] [pat_1] [pat_2] [insertSize] [insertVar] [ref]" << endl << endl
         << "Reformat paired-end alignments into SAM output." << endl
         << "\t[al_1]\t\tOutput from gcsa_alignment (paired-ends in separate files)" << endl
         << "\t[al_2]" << endl
         << "\t[pat_1]\t\tPaired-end reads in separate text files" << endl
         << "\t[pat_2]" << endl
         << "\t[insertSize]\tEstimated insert size" << endl
         << "\t[insertVar]\tEstimated insert size variance" << endl
         << "\t[ref]\t\tThe reference sequence (as a FASTA file)" << endl << endl
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

void myassert(char *p, char *t)
{
    if (strcmp(p,t) == 0)
        return;
    revstr(p);
    complement(p);
    if (strcmp(p,t) != 0)
    {
        cerr << "myassert() failed p and t do not match" << endl;
        cerr << p << endl;
        cerr << t << endl;
        exit(2);
    }
    
    // Revert the string
    revstr(p);
    complement(p);
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

pair<unsigned,pair<int,unsigned> > rescue(string ref, char *pattern, vector<int> const &pos, int insertSize, int insertVar)
{
    unsigned minedit = strlen(pattern);
    unsigned mini = 0;
    int minpos = 0;
    MyersEditDistance ed((uchar *)pattern, strlen(pattern), strlen(pattern));
    for(unsigned i = 0; i < pos.size(); ++i)
    {
        // Check both sides of each position
        int p = pos[i];
        string tmp = ref.substr(p + insertSize - insertVar, 2*insertVar + strlen(pattern));
        pair<unsigned,unsigned> r = ed.prefixDist((uchar *)tmp.c_str(), tmp.size());
        if (minedit > r.first )
        {
            minedit = r.first;
            minpos = p + insertSize - insertVar + r.second;
            mini = i;
        }
        
        if (p > insertSize + insertVar)
        { 
            tmp = ref.substr(p - insertSize - insertVar, 2*insertVar + strlen(pattern));
            //cerr << tmp << endl;
            //cerr << pattern << endl;
            r = ed.prefixDist((uchar *)tmp.c_str(), tmp.size());
            if (minedit > r.first )
            {
                minedit = r.first;
                minpos = p - insertSize - insertVar + r.second;
                mini = i;
            }
        }
    }

#ifdef DEBUG
    cerr << "rescuing " << pattern << " with ed = " << minedit << " at " << minpos << endl;
#endif
    return make_pair(minedit, make_pair(minpos, mini));
}

int main(int argc, char **argv) 
{
    if (argc != 8)
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
    FILE *iss = fopen(argv[2], "r");
    if (iss == 0)
    {
        cerr << "error: al_2 file " << argv[2] << " not found" << endl;
        return 1;
    }
    FILE *patternf = fopen(argv[3], "r");
    if (patternf == 0)
    {
        cerr << "error: pat_1 file " << argv[3] << " not found" << endl;
        return 1;
    }
    FILE *patterns = fopen(argv[4], "r");
    if (patterns == 0)
    {
        cerr << "error: pat_2 file " << argv[4] << " not found" << endl;
        return 1;
    }
    int insertSize = atoi(argv[5]);
    int insertVar = atoi(argv[6]);
    FILE *reff = fopen(argv[7], "r");
    if (reff == 0)
    {
        cerr << "error: reference file " << argv[7] << " not found" << endl;
        return 1;
    }
    string ref = readRef(reff);
    //cerr << "ref: \"" << ref << "\"" << endl;;
    fclose(reff);

    int f = -1, s = -1; // current read id
    char readf[L], reads[L]; // current read pattern
    vector<int> posf, poss; // vector of current read positions
    char cachef[L], caches[L]; // first row of the next pattern
    fgets(cachef, L, isf);
    readNextRead(isf, f, readf, posf, cachef);
    fgets(caches, L, iss);
    readNextRead(iss, s, reads, poss, caches);

    char patf[L], pats[L];
    int nextOutput = 1;
    unsigned rescued = 0, failed = 0, rescuedhit = 0, bothhit = 0;

    while (!feof(isf) || !feof(iss))
    {        
        // Output unmapped pairs
        while (nextOutput < s && nextOutput < f)
        {
            cout << nextOutput <<  "\t77\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t0\t0\t*\t*\n"; 
            cout << nextOutput << "\t141\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t0\t0\t*\t*\n"; 
            fgets(patf, L, patternf); // Skip in pattern file
            fgets(pats, L, patterns);
            ++nextOutput;
        }

        fgets(patf, L, patternf); // Keep pattern files in sync
        fgets(pats, L, patterns);
        patf[strlen(patf)-1]=0; 
        pats[strlen(pats)-1]=0;

        if (f < s)
        {
            // First of pair was mapped (mate unmapped)
#ifdef DEBUG
            cerr << "id = " << f << ", read = " << readf << ", pos =";
            for (vector<int>::iterator it = posf.begin(); it != posf.end(); ++it)
                cerr << ", " << *it;
            cerr << endl;
#endif
            myassert(readf, patf);

            // Mate rescue for second in pair
            pair<unsigned, pair<int,unsigned> > forw = rescue(ref, pats, posf, insertSize, insertVar);
            revstr(pats); complement(pats);
            pair<unsigned, pair<int,unsigned> > rc = rescue(ref, pats, posf, insertSize, insertVar);
            revstr(pats); complement(pats);
            if (forw.first > rc.first)
            {
                forw.first = rc.first;
                forw.second = rc.second;
            }
            if (forw.first > 12)
                forw.second.first = 0;
            pair<int,unsigned> p = forw.second;

            if (p.first > 0)
            {
                // Successful rescue
                ++rescued;
                cout << f << "\t67\tgi|224589812|ref|nc_000020.10|\t" << posf[p.second] << "\t255\t*\t*\t" << p.first << "\t0\t*\t*\n"; 
                cout << f << "\t131\tgi|224589812|ref|nc_000020.10|\t" << p.first << "\t255\t*\t*\t" << posf[p.second] << "\t0\t*\t*\n"; 
            }
            else
            {
                // Still unmapped
                ++failed;
                cout << f << "\t73\tgi|224589812|ref|nc_000020.10|\t" << posf[0] << "\t255\t*\t*\t0\t0\t*\t*\n"; 
                cout << f << "\t133\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t" << posf[0] << "\t0\t*\t*\n"; 
            }

            readNextRead(isf, f, readf, posf, cachef);
        }
        else if (f > s)
        {
            // Second of pair was mapped (mate unmapped)
#ifdef DEBUG
            cerr << "id = " << s << ", read = " << reads << ", pos =";
            for (vector<int>::iterator it = poss.begin(); it != poss.end(); ++it)
                cerr << ", " << *it;
            cerr << endl;
#endif
            myassert(reads, pats);

            // Mate rescue for first in pair
            pair<unsigned, pair<int,unsigned> > forw = rescue(ref, patf, poss, insertSize, insertVar);
            revstr(patf); complement(patf);
            pair<unsigned, pair<int,unsigned> > rc = rescue(ref, patf, poss, insertSize, insertVar);
            revstr(patf); complement(patf);
            if (forw.first > rc.first)
            {
                forw.first = rc.first;
                forw.second = rc.second;
            }
            if (forw.first > 12)
                forw.second.first = 0;
            pair<int,unsigned> p = forw.second;

            if (p.first > 0)
            {
                // Successful rescue
                ++rescued;
                cout << s << "\t67\tgi|224589812|ref|nc_000020.10|\t" << p.first << "\t255\t*\t*\t" << poss[p.second] << "\t0\t*\t*\n"; 
                cout << s << "\t131\tgi|224589812|ref|nc_000020.10|\t" << poss[p.second] << "\t255\t*\t*\t" << p.first << "\t0\t*\t*\n"; 
            }
            else
            {
                // Still unmapped
                ++failed;
                cout << s << "\t69\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t" << poss[0] << "\t0\t*\t*\n"; 
                cout << s << "\t137\tgi|224589812|ref|nc_000020.10|\t" << poss[0] << "\t255\t*\t*\t0\t0\t*\t*\n"; 
            }

            readNextRead(iss, s, reads, poss, caches);
        }
        else
        {
            // Both mapped
            assert(f == s);
            myassert(readf, patf);
            myassert(reads, pats);
            bothhit++;

            // Choose the matepair to output: take the one having distance closest to insertSize
            int mindist = 100*insertSize;
            unsigned minf = 0, mins = 0;
            for (unsigned i = 0; i < posf.size(); ++i)
                for (unsigned j = 0; j < poss.size(); ++j)
                    if (mindist > abs(insertSize-abs(posf[i]-poss[j])))
                    {
                        mindist = abs(insertSize-abs(posf[i]-poss[j]));
                        minf = i;
                        mins = j;
                    }

            {
                // Mate rescue for second in pair
                pair<unsigned, pair<int,unsigned> > forw = rescue(ref, pats, posf, insertSize, insertVar);
                revstr(pats); complement(pats);
                pair<unsigned, pair<int,unsigned> > rc = rescue(ref, pats, posf, insertSize, insertVar);
                revstr(pats); complement(pats);
                if (forw.first > rc.first)
                {
                    forw.first = rc.first;
                    forw.second = rc.second;
                }
                if (forw.first > 8)
                    forw.second.first = 0;
                pair<int,unsigned> p = forw.second;
                
                if (p.first > 0 && mindist > abs(insertSize-abs(posf[p.second]-p.first)))
                {
                    // Successful rescue of second in pair
                    mindist = abs(insertSize-abs(posf[p.second]-p.first));
                    minf = p.second;
                    mins = 0;
                    poss[mins] = p.first;
                    rescuedhit++;
                }
            }
            {
                // Mate rescue for first in pair
                pair<unsigned, pair<int,unsigned> > forw = rescue(ref, patf, poss, insertSize, insertVar);
                revstr(patf); complement(patf);
                pair<unsigned, pair<int,unsigned> > rc = rescue(ref, patf, poss, insertSize, insertVar);
                revstr(patf); complement(patf);
                if (forw.first > rc.first)
                {
                    forw.first = rc.first;
                    forw.second = rc.second;
                }
                if (forw.first > 8)
                    forw.second.first = 0;
                pair<int,unsigned> p = forw.second;
                
                if (p.first > 0 && abs(insertSize-abs(poss[p.second]-p.first)))
                {
                    // Successful rescue
                    mindist = abs(insertSize-abs(poss[p.second]-p.first));
                    mins = p.second; 
                    minf = 0;
                    posf[minf]=p.first;
                    rescuedhit++;
                }
            }
            
            cout << f << "\t67\tgi|224589812|ref|nc_000020.10|\t" << posf[minf] << "\t255\t*\t*\t" << poss[mins] << "\t0\t*\t*\n"; 
            cout << s << "\t131\tgi|224589812|ref|nc_000020.10|\t" << poss[mins] << "\t255\t*\t*\t" << posf[minf] << "\t0\t*\t*\n"; 

            readNextRead(isf, f, readf, posf, cachef);
            readNextRead(iss, s, reads, poss, caches);
        }
        nextOutput ++;
        if (nextOutput % 100000 == 0)
            cerr << "nextOutput = " << nextOutput << ", " << rescued << " read pairs rescued, " << rescuedhit << " hits rescued, with <8e (" << 100*rescued/(rescued+failed) << "%, " << 100*rescuedhit/bothhit << "%)." << endl;
    }

    // Flush remaining unmapped pairs
    while (fgets(patf,L,patternf))
    {
        cout << nextOutput <<  "\t77\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t0\t0\t*\t*\n"; 
        cout << nextOutput << "\t141\tgi|224589812|ref|nc_000020.10|\t0\t255\t*\t*\t0\t0\t*\t*\n"; 
        ++nextOutput;
    }

    cerr << "All input files exausted. " << rescued << " read pairs rescued." << endl;
    return 0;
}
