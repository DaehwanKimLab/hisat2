// Myers edit distance originally implemented by Veli M‰kinen
#ifndef _Tools_h_
#define _Tools_h_

#include <utility>
#include <climits>

typedef unsigned long ulong;
typedef unsigned char uchar;
#define W (CHAR_BIT*sizeof(unsigned long))

class EditDistance 
{
public:
    uchar const *pat;
    unsigned  patlen;
    unsigned maxlen;
    unsigned **dp;
    unsigned p;

    bool end()
    {
        if (p < maxlen)
            return false;
        return true;
    }

    EditDistance(uchar const *pattern, unsigned len, unsigned mlen) {
        pat = pattern;
        patlen = len;
        maxlen = mlen;
        p = 0;
        dp = new unsigned*[patlen + 1];
        for (unsigned i = 0; i <= patlen; i++) {
            dp[i] = new unsigned[maxlen+1];
            dp[i][0] = i;
        }
    }

    ~EditDistance() {
        for (unsigned i = 0; i <= patlen; i++)
            delete [] dp[i];
        delete [] dp;
    }

    std::pair<unsigned, unsigned> prefixDist(uchar const *cand, unsigned candlen) {
        unsigned mindist = patlen+candlen;
        unsigned minlen = patlen+candlen;
        unsigned *col = new unsigned[patlen+1];
        unsigned *ncol = new unsigned[patlen+1];
        for (unsigned i = 0; i <= patlen; i++)
            col[i] = i;
        for (unsigned j = 1; j <= candlen; j++) {
            ncol[0] = 0;
            for (unsigned i = 1; i <= patlen; i++) {
                unsigned match = (cand[j-1] == pat[i-1])? col[i - 1] : (col[i - 1] + 1);
                ncol[i] = std::min(col[i] + 1, ncol[i - 1] + 1);
                ncol[i] = std::min(ncol[i],match);
            }
            if (ncol[patlen]<mindist)
            {
                mindist = ncol[patlen];
                minlen = j;
            }
            for (unsigned i = 0; i <= patlen; i++)
                col[i] = ncol[i];
        }
        delete [] col;
        delete [] ncol;
        return std::make_pair(mindist,minlen);
    }


    unsigned pushChar(uchar c) {
        unsigned mindist = patlen;
        p++;
        if (p>maxlen)
            std::cout << "out of bounds in DP computation, p = " << p << ", max = " << maxlen << std::endl;
        dp[0][p] = p;
        for (unsigned i = 1; i <= patlen; i++) {
            unsigned match = (c == pat[i-1])? dp[i - 1][p-1] : (dp[i - 1][p-1] + 1);
            dp[i][p] = std::min(dp[i][p-1] + 1, dp[i - 1][p] + 1);
            dp[i][p] = std::min(dp[i][p],match);
            if (dp[i][p]<mindist)
                mindist = dp[i][p];
        }
        return mindist;
    }

    void popChar() {
        --p;
        return;
    }

    unsigned dist() {
        return dp[patlen][p];
    }

    unsigned getSuffixLength()
    {
        unsigned mindist = dp[patlen][p];
        unsigned minlen = patlen;
        for (unsigned i = 0; i < patlen; ++i)
            if (dp[i][p] < mindist)
            {
                mindist = dp[i][p];
                minlen = i;
            }
        return minlen;
    }

};


class MyersEditDistance
{
public:
    ulong **Peq;
    bool S[256];
    unsigned k;
    unsigned m;
    unsigned blockCount;
    unsigned block;
    unsigned lastBlockLength;
    ulong *Pv;
    ulong *Mv;
    ulong *Score;
    ulong twopowmpos;
    ulong twopowwpos;
 
    void MPreComp(uchar const *P, unsigned m)
    {
        this->m = m;
        ulong i;
        ulong j;
        Peq = new ulong *[(m-1) / W + 1]; // (unsigned **)malloc((((m-1) / w)+1)*sizeof(unsigned*));
        ulong l;
        for (l=0;l <= ((m-1) / W); l++)
            Peq[l] = new ulong[256]; //(unsigned *)malloc(256*sizeof(unsigned));
        for(i = 0;i< 256;i++) {
            S[i] = false;
            for (l=0;l <= ((m-1) / W); l++)
                Peq[l][i] = 0;
        }
        for (l=0;l <= ((m-1) / W); l++)
            for (j = 1; j <= W; j++) {
                if (W*l+j > m)
                    break;
                Peq[l][P[W*l+j-1]] = Peq[l][P[W*l+j-1]] | (1lu << (j-1));
                if (j+W*l <= k+1)
                    S[(unsigned)P[W*l+j-1]] = true;
            }
    }

    MyersEditDistance(uchar const *P, unsigned m, unsigned maxk) 
    : Peq(0), k(maxk)
    {
        MPreComp(P, m);
        blockCount= ((m-1) / W);
        block = (k-1) / W;
        lastBlockLength = m % W;
        if (lastBlockLength == 0)
            lastBlockLength = W;
        Pv = new ulong[blockCount+1];  // (unsigned *)malloc((blockCount+1)*sizeof(unsigned));
        Mv = new ulong[blockCount+1];  // (unsigned *)malloc((blockCount+1)*sizeof(unsigned));
        Score = new ulong[blockCount+1];//(unsigned *)malloc((blockCount+1)*sizeof(unsigned));
        twopowmpos = 1lu << (lastBlockLength-1);
        twopowwpos = 1lu << (W-1);
        ulong l;
        for (l=0; l < blockCount; l++) {
            Mv[l] = 0lu;
            Score[l] = (l+1lu)*W;
            Pv[l] = ~0lu;
        }
        Mv[blockCount] = 0;
        Score[blockCount] = m;
        Pv[blockCount] = (~0lu) >> ((blockCount+1lu)*W - m);
    }

    ~MyersEditDistance()
    {
        delete [] Pv; Pv = 0;
        delete [] Mv; Mv = 0;
        delete [] Score; Score = 0;
        unsigned l;
        for (l=0;l <= ((m-1) / W); l++)
        {
            delete [] Peq[l];
            Peq[l] = 0;
        }
        delete [] Peq; Peq = 0;
    }

    // FIXME if m < W use the simpler computation
    std::pair<unsigned, unsigned> prefixDist(uchar const *text, unsigned n)
    {
        unsigned mindist = m+n;
        unsigned minlen = m+n;
        //            int count=0;
        //int blockAverage=0;
        ulong Eq, Xv, Xh, Ph, Mh;
        ulong Phendbit;
        ulong Mhendbit;
        ulong temp;
        ulong j;
        ulong l;
        ulong overk=1;
        block = (k-1) / W;
        for(j = 1; j <= n; j++) {
            Phendbit = 0; //1lu; //0;
            Mhendbit = 0;
            for (l=0; l <= blockCount; l++) {
                Eq = Peq[l][(unsigned char)text[j-1]];
                Xv = Eq | Mv[l];
                Xh = ((((Eq | Mhendbit)& Pv[l]) + Pv[l]) ^ Pv[l]) | (Eq | Mhendbit);
                Ph = Mv[l] | ~ (Xh | Pv[l]);
                Mh = Pv[l] & Xh;
                temp = l < blockCount ? twopowwpos : twopowmpos;
                if (Ph & temp)
                    Score[l] += 1;
                else if (Mh & temp)
                    Score[l] -= 1;
                temp = (Ph >> (W-1));
                Ph <<= 1;
                Ph |= Phendbit;
                Phendbit = temp;
                temp = (Mh >> (W-1));
                Mh <<= 1;
                Mh |= Mhendbit;
                Mhendbit = temp;
                Pv[l] = (Mh) | ~ (Xv | Ph);
                Mv[l] = Ph & Xv;
                if (block == l+1 &&
                    ((block == blockCount &&
                      (Score[l] > k+lastBlockLength ||
                       (Score[l] > k && overk != 0 && j - overk >= lastBlockLength))) ||
                     (block != blockCount &&
                      (Score[l] > k+W ||
                       (Score[l] > k && overk != 0 && j - overk >= W)))))  {
                    // lohkon Score kasvanut niin isoksi, ettei seuraavaa kannata laskea
                    overk = 0;
                    block = l;
                    break;
                }
                else if (block == l+1 && overk == 0 && Score[l] > k)
                    // Talletetaan ensimm‰inen diagonaali jolla Score > k. N‰in tiedet‰‰n
                    // jatkossa koska seuraavan lohkon laskeminen on turhaa.
                    overk = j;
                else if (block == l+1 && Score[l] <= k)
                    // Diagonaalilla Score <= k ==> Nollataan muuttuja overk.
                    overk = 0;
                else if (block == l && block != blockCount && Score[l] <= k) {
                    // Score <= k ==> jatketaan seuraavaan lohkoon. Koska seuraavaa lohkoa
                    // ei ole laskettu edellisell‰ sarakkeella, pit‰‰ sen arvot p‰ivitt‰‰.
                    Score[l+1] = k+1lu;
                    Pv[l+1] = 0lu;
                    Mv[l+1] = 0lu;
                    overk = 0lu;
                    block = l+1;
                }
                else if (block == l)
                    // Seuraavaan lohkoon ei kannata siirty‰, kun Score > k.
                    break;
            }
//                blockAverage += block+1;
            // if (block == blockCount && Score[blockCount] <= k) {
            if (block == blockCount && Score[blockCount] < mindist)
            {
                mindist = Score[blockCount];
                minlen = j;
                
                //        printf("Match at %d\n",j);
                //    count++;
            }
        }
//            return count;
        return std::make_pair(mindist,minlen);
        /*ShowMessage("Esiintymi‰ lˆytyi " + IntToStr(laskuri) + "kpl. " +
          "Aikaa meni " + IntToStr(msec) + "msec");
          ShowMessage("Keskim‰‰r‰inen lohkojen m‰‰r‰ " +
          FormatFloat("0.00",(float)blockAverage / n));*/
    }

};

#endif
