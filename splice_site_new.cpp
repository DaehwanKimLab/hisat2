/*
 * Copyright 2013, Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "edit.h"
#include "splice_site.h"
#include "aligner_report.h"
#include "aligner_result.h"

#if defined(NEW_PROB_MODEL)

#include "splice_site_mem.h"

#else

float donor_prob[4][donor_len] = {
    {0.340f, 0.604f, 0.092f, 0.001f, 0.001f, 0.526f, 0.713f, 0.071f, 0.160f},
    {0.363f, 0.129f, 0.033f, 0.001f, 0.001f, 0.028f, 0.076f, 0.055f, 0.165f},
    {0.183f, 0.125f, 0.803f, 1.000f, 0.001f, 0.419f, 0.118f, 0.814f, 0.209f},
    {0.114f, 0.142f, 0.073f, 0.001f, 1.000f, 0.025f, 0.093f, 0.059f, 0.462f}
};

float acceptor_prob[4][acceptor_len] = {
    {0.090f, 0.084f, 0.075f, 0.068f, 0.076f, 0.080f, 0.097f, 0.092f, 0.076f, 0.078f, 0.237f, 0.042f, 1.000f, 0.001f, 0.239f},
    {0.310f, 0.310f, 0.307f, 0.293f, 0.326f, 0.330f, 0.373f, 0.385f, 0.410f, 0.352f, 0.309f, 0.708f, 0.001f, 0.001f, 0.138f},
    {0.125f, 0.115f, 0.106f, 0.104f, 0.110f, 0.113f, 0.113f, 0.085f, 0.066f, 0.064f, 0.212f, 0.003f, 0.001f, 1.000f, 0.520f},
    {0.463f, 0.440f, 0.470f, 0.494f, 0.471f, 0.463f, 0.408f, 0.429f, 0.445f, 0.504f, 0.240f, 0.246f, 0.001f, 0.001f, 0.104f}
};

float donor_prob_sum[1 << (donor_len << 1)];
float acceptor_prob_sum1[1 << (acceptor_len1 << 1)];
float acceptor_prob_sum2[1 << (acceptor_len2 << 1)];

#endif

void init_junction_prob()
{
#if !defined(NEW_PROB_MODEL)
    for(size_t i = 0; i < donor_len; i++) {
        ASSERT_ONLY(float sum = 0.0f);
        for(size_t j = 0; j < 4; j++) {
            float prob = donor_prob[j][i];
            assert_gt(prob, 0.0f);
            ASSERT_ONLY(sum += prob);
            donor_prob[j][i] = log(prob / background_prob[j]);
        }
        assert_range(0.9f, 1.05f, sum);
    }
    for(size_t i = 0; i < acceptor_len; i++) {
        ASSERT_ONLY(float sum = 0.0f);
        for(size_t j = 0; j < 4; j++) {
            float prob = acceptor_prob[j][i];
            assert_gt(prob, 0.0f);
            ASSERT_ONLY(sum += prob);
            acceptor_prob[j][i] = log(prob / background_prob[j]);
        }
        assert_range(0.9f, 1.05f, sum);
    }
    
    const size_t donor_elms = 1 << (donor_len << 1);
    for(size_t i = 0; i < donor_elms; i++) {
        float sum = 0.0f;
        for(size_t j = 0; j < donor_len; j++) {
            int base = (i >> (j << 1)) & 0x3;
            sum += donor_prob[base][donor_len - j - 1];
        }
        donor_prob_sum[i] = exp(-sum);
    }
    
    const size_t acceptor_elms1 = 1 << (acceptor_len1 << 1);
    for(size_t i = 0; i < acceptor_elms1; i++) {
        float sum = 0.0f;
        for(size_t j = 0; j < acceptor_len1; j++) {
            int base = (i >> (j << 1)) & 0x3;
            sum += acceptor_prob[base][acceptor_len1 - j - 1];
        }
        acceptor_prob_sum1[i] = exp(-sum);
    }
    
    const size_t acceptor_elms2 = 1 << (acceptor_len2 << 1);
    for(size_t i = 0; i < acceptor_elms2; i++) {
        float sum = 0.0f;
        for(size_t j = 0; j < acceptor_len2; j++) {
            int base = (i >> (j << 1)) & 0x3;
            sum += acceptor_prob[base][acceptor_len - j - 1];
        }
        acceptor_prob_sum2[i] = exp(-sum);
    }
#endif
}

ostream& operator<<(ostream& out, const SpliceSite& s)
{
	out << s.ref() << "\t"
        << s.left() << "\t"
        << s.right() << "\t"
        << (s.fw() ? "+" : "-") << endl;
	return out;
}

#if defined(SPLIT_DB)

static const uint64_t ref_segment_size = 1 << 24;

SpliceSiteDB::SpliceSiteDB(
                           const BitPairReference& refs,
                           const EList<string>& refnames,
                           bool threadSafe,
                           bool write,
                           bool read) :
_numRefs(refs.numRefs()),
_refnames(refnames),
_write(write),
_read(read),
_threadSafe(threadSafe),
_empty(true)
{
    assert_gt(_numRefs, 0);
    assert_eq(_numRefs, _refnames.size());
    for(uint64_t i = 0; i < _numRefs; i++) {
      uint64_t reflen = refs.approxLen(i);
      assert_gt(reflen, 0);
      uint64_t numsegs = (reflen + ref_segment_size - 1) / ref_segment_size;
      assert_gt(numsegs, 0);
      
      _fwIndex.expand();
      _bwIndex.expand();
      for(uint64_t j = 0; j < numsegs; j++) {
	_fwIndex.back().push_back(new RedBlack<SpliceSitePos, uint32_t>(16 << 10, CA_CAT));
	_bwIndex.back().push_back(new RedBlack<SpliceSitePos, uint32_t>(16 << 10, CA_CAT));
      }

      _pool.expand();
      _pool.back().resize(numsegs);
      _spliceSites.expand();
      _spliceSites.back().resize(numsegs);

      _mutex.expand();
      for(uint64_t j = 0; j < numsegs; j++) {
	_mutex.back().push_back(MUTEX_T());
      }
    }
    
    donorstr.resize(donor_exonic_len + donor_intronic_len);
    acceptorstr.resize(acceptor_intronic_len + acceptor_exonic_len);
}

SpliceSiteDB::~SpliceSiteDB() {
    assert_eq(_fwIndex.size(), _bwIndex.size());
    assert_eq(_fwIndex.size(), _pool.size());
    for(uint64_t i = 0; i < _numRefs; i++) {
      assert_eq(_fwIndex[i].size(), _bwIndex[i].size());
      for(uint64_t j = 0; j < _fwIndex[i].size(); j++) {
        delete _fwIndex[i][j];
        delete _bwIndex[i][j];
      }

      _fwIndex[i].clear();
      _bwIndex[i].clear();
        
        ELList<Pool*>& pool = _pool[i];
        for(size_t j = 0; j < pool.size(); j++) {
	  for(size_t k = 0; k < pool[j].size(); k++) {
            delete pool[j][k];
	  }
        }
    }
}

#if 0
size_t SpliceSiteDB::size(uint64_t ref) const {
    if(!_read) return 0;

    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    assert_lt(ref, _fwIndex.size());
    assert_eq(_fwIndex.size(), _bwIndex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    return _fwIndex.size();
}

bool SpliceSiteDB::empty(uint64_t ref) const {
    return size(ref) == 0;
}
#endif

bool SpliceSiteDB::addSpliceSite(
                                 const Read& rd,
                                 const AlnRes& rs,
                                 uint32_t minAnchorLen)
{
    if(!_write) return false;
    if(rs.trimmed5p(true) + rs.trimmed3p(true) > 0) return false;
    
    _empty = false;
    
    Coord coord = rs.refcoord();
    uint64_t ref = coord.ref();
    assert_lt(ref, _numRefs);
    const EList<Edit>& edits = rs.ned();
    if(!coord.orient()) {
        Edit::invertPoss(const_cast<EList<Edit>&>(edits), rd.length(), false);
	}
    SpliceSitePos ssp;
    uint32_t refoff = coord.off();
    uint32_t leftAnchorLen = 0, rightAnchorLen = 0;
    size_t eidx = 0;
    size_t last_eidx = 0;
    for(size_t i = 0; i < rd.length(); i++, refoff++) {
        while(eidx < edits.size() && edits[eidx].pos == i) {
			if(edits[eidx].isReadGap()) {
                refoff++;
			} else if(edits[eidx].isRefGap()) {
                assert_gt(refoff, 0);
                refoff--;
            }
            if(edits[eidx].isSpliced()) {
                assert_gt(refoff, 0);
                if(ssp.inited()) {
                    assert(edits[last_eidx].isSpliced());
                    assert_lt(edits[last_eidx].pos, edits[eidx].pos);
                    rightAnchorLen = edits[eidx].pos - edits[last_eidx].pos;
                    if(leftAnchorLen >= minAnchorLen && rightAnchorLen >= minAnchorLen) {
                        bool added = false;
                        uint64_t ref_seg = ssp.left() / ref_segment_size;
                        uint64_t ref_seg2 = ssp.right() / ref_segment_size;
                        assert_leq(ref_seg, ref_seg2);
                        assert_lt(ref, _mutex.size());
                        assert_lt(ref_seg, _mutex[ref].size());
                        assert_lt(ref_seg2, _mutex[ref].size());
                        ThreadSafe t(&_mutex[ref][ref_seg], _threadSafe && _write);
                        ThreadSafe t2(&_mutex[ref][ref_seg2], _threadSafe && _write && ref_seg != ref_seg2);
                        assert_lt(ref, _fwIndex.size());
                        assert_lt(ref_seg, _fwIndex[ref].size());
                        assert(_fwIndex[ref][ref_seg] != NULL);
                        Node *cur = _fwIndex[ref][ref_seg]->add(pool(ref, ref_seg), ssp, &added);
                        if(added) {
                            assert_lt(ref, _spliceSites.size());
                            assert_lt(ref_seg, _spliceSites[ref].size());
                            _spliceSites[ref][ref_seg].expand();
                            _spliceSites[ref][ref_seg].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                            _spliceSites[ref][ref_seg].back()._readid = rd.rdid;
                            assert(cur != NULL);
                            cur->payload = _spliceSites[ref][ref_seg].size() - 1;
                            
                            SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                            assert_lt(ref, _bwIndex.size());
                            assert_lt(ref_seg2, _bwIndex[ref].size());
                            assert(_bwIndex[ref][ref_seg2] != NULL);
                            cur = _bwIndex[ref][ref_seg2]->add(pool(ref, ref_seg2), rssp, &added);
                            assert(added);
                            assert(cur != NULL);
                            if(ref_seg != ref_seg2) {
                                _spliceSites[ref][ref_seg2].expand();
                                _spliceSites[ref][ref_seg2].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                                _spliceSites[ref][ref_seg2].back()._readid = rd.rdid;
                            }
                            cur->payload = _spliceSites[ref][ref_seg2].size() - 1;
                            
                        } else {
                            assert(cur != NULL);
                            assert_lt(ref, _spliceSites.size());
                            assert_lt(ref_seg, _spliceSites[ref].size());
                            assert_lt(cur->payload, _spliceSites[ref][ref_seg].size());
                            if(rd.rdid < _spliceSites[ref][ref_seg][cur->payload]._readid) {
                                _spliceSites[ref][ref_seg][cur->payload]._readid = rd.rdid;
                            }
                            if(ref_seg != ref_seg2) {
                                SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                                cur = _bwIndex[ref][ref_seg2]->add(pool(ref, ref_seg2), rssp, &added);
                                assert(cur != NULL);
                                assert(!added);
                                if(rd.rdid < _spliceSites[ref][ref_seg2][cur->payload]._readid) {
                                    _spliceSites[ref][ref_seg2][cur->payload]._readid = rd.rdid;
                                }
                            }
                        }
                    }
                    leftAnchorLen = rightAnchorLen;
                    rightAnchorLen = 0;
                } else {
                    leftAnchorLen = edits[eidx].pos;
                }
                bool fw = (edits[eidx].splDir == EDIT_SPL_FW || edits[eidx].splDir == EDIT_SPL_UNKNOWN);
                bool canonical = (edits[eidx].splDir != EDIT_SPL_UNKNOWN);
                ssp.init(coord.ref(), refoff - 1, refoff + edits[eidx].splLen, fw, canonical);
                refoff += edits[eidx].splLen;
                last_eidx = eidx;
            }
			eidx++;
		}
	}
    if(ssp.inited()) {
        assert(edits[last_eidx].isSpliced());
        assert_lt(edits[last_eidx].pos, rd.length());
        rightAnchorLen = rd.length() - edits[last_eidx].pos;
        if(leftAnchorLen >= minAnchorLen && rightAnchorLen >= minAnchorLen) {
            bool added = false;
            uint64_t ref_seg = ssp.left() / ref_segment_size;
            uint64_t ref_seg2 = ssp.right() / ref_segment_size;
            assert_leq(ref_seg, ref_seg2);
            assert_lt(ref, _mutex.size());
            assert_lt(ref_seg, _mutex[ref].size());
            assert_lt(ref_seg2, _mutex[ref].size());
            ThreadSafe t(&_mutex[ref][ref_seg], _threadSafe && _write);
            ThreadSafe t2(&_mutex[ref][ref_seg2], _threadSafe && _write && ref_seg != ref_seg2);
            assert_lt(ref, _fwIndex.size());
            assert_lt(ref_seg, _fwIndex[ref].size());
            assert(_fwIndex[ref][ref_seg] != NULL);
            Node *cur = _fwIndex[ref][ref_seg]->add(pool(ref, ref_seg), ssp, &added);
            if(added) {
                assert_lt(ref, _spliceSites.size());
                assert_lt(ref_seg, _spliceSites[ref].size());
                _spliceSites[ref][ref_seg].expand();
                _spliceSites[ref][ref_seg].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                _spliceSites[ref][ref_seg].back()._readid = rd.rdid;
                assert(cur != NULL);
                cur->payload = _spliceSites[ref][ref_seg].size() - 1;
                
                SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                assert_lt(ref, _bwIndex.size());
                assert_lt(ref_seg2, _bwIndex[ref].size());
                assert(_bwIndex[ref][ref_seg2] != NULL);
                cur = _bwIndex[ref][ref_seg2]->add(pool(ref, ref_seg2), rssp, &added);
                assert(added);
                assert(cur != NULL);
                if(ref_seg != ref_seg2) {
                    _spliceSites[ref][ref_seg2].expand();
                    _spliceSites[ref][ref_seg2].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                    _spliceSites[ref][ref_seg2].back()._readid = rd.rdid;
                }
                cur->payload = _spliceSites[ref][ref_seg2].size() - 1;
                
            } else {
                assert(cur != NULL);
                assert_lt(ref, _spliceSites.size());
                assert_lt(ref_seg, _spliceSites[ref].size());
                assert_lt(cur->payload, _spliceSites[ref][ref_seg].size());
                if(rd.rdid < _spliceSites[ref][ref_seg][cur->payload]._readid) {
                    _spliceSites[ref][ref_seg][cur->payload]._readid = rd.rdid;
                }
                if(ref_seg != ref_seg2) {
                    SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                    cur = _bwIndex[ref][ref_seg2]->add(pool(ref, ref_seg2), rssp, &added);
                    assert(cur != NULL);
                    assert(!added);
                    if(rd.rdid < _spliceSites[ref][ref_seg2][cur->payload]._readid) {
                        _spliceSites[ref][ref_seg2][cur->payload]._readid = rd.rdid;
                    }
                }
            }
        }
    }
    if(!coord.orient()) {
		Edit::invertPoss(const_cast<EList<Edit>&>(edits), rd.length(), false);
	}

    return true;
}

bool SpliceSiteDB::getSpliceSite(SpliceSite& ss) const
{
    if(!_read) return false;

    uint64_t ref = ss.ref();
    uint64_t ref_seg = ss.left() / ref_segment_size;
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    assert_lt(ref_seg, _mutex[ref].size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref][ref_seg]), _threadSafe && _write);
    
    assert_lt(ref, _fwIndex.size());
    assert_lt(ref_seg, _fwIndex[ref].size());
    assert(_fwIndex[ref][ref_seg] != NULL);
    const Node *cur = _fwIndex[ref][ref_seg]->lookup(ss);
    if(cur == NULL) return false;
    assert(cur != NULL);
    assert_lt(ref, _spliceSites.size());
    assert_lt(ref_seg, _spliceSites[ref].size());
    ss = _spliceSites[ref][ref_seg][cur->payload];
    return true;
}

void SpliceSiteDB::getLeftSpliceSites(uint32_t ref, uint32_t left, uint32_t range, EList<SpliceSite>& spliceSites) const
{
    if(!_read) return;
    spliceSites.clear();
    assert_lt(ref, _numRefs);
    assert_gt(range, 0);
    assert_geq(left + 1, range);
    uint32_t begin = left + 1 - range, end = left;
    for(uint32_t i = begin / ref_segment_size; i <= end / ref_segment_size; i++) {
        assert_lt(ref, _mutex.size());
        assert_lt(i, _mutex[ref].size());
        ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref][i]), _threadSafe && _write);
        assert_lt(ref, _bwIndex.size());
        assert_lt(i, _bwIndex[ref].size());
        assert(_bwIndex[ref][i] != NULL);
        const Node *cur = _bwIndex[ref][i]->root();
        if(cur != NULL) getSpliceSites_recur(cur, _spliceSites[ref][i], ref, begin, end, spliceSites);
    }
}

void SpliceSiteDB::getRightSpliceSites(uint32_t ref, uint32_t right, uint32_t range, EList<SpliceSite>& spliceSites) const
{
    if(!_read) return;
    spliceSites.clear();
    assert_lt(ref, _numRefs);
    assert_gt(range, 0);
    assert_gt(right + range, range);
    uint32_t begin = right, end = right + range - 1;
    for(uint32_t i = begin / ref_segment_size; i <= end / ref_segment_size; i++) {
        assert_lt(ref, _mutex.size());
        assert_lt(i, _mutex[ref].size());
        ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref][i]), _threadSafe && _write);
        assert_lt(ref, _fwIndex.size());
        assert_lt(i, _fwIndex[ref].size());
        assert(_fwIndex[ref][i] != NULL);
        const Node *cur = _fwIndex[ref][i]->root();
        if(cur != NULL) getSpliceSites_recur(cur, _spliceSites[ref][i], ref, begin, end, spliceSites);
    }
    
}

void SpliceSiteDB::getSpliceSites_recur(
                                        const RedBlackNode<SpliceSitePos, uint32_t> *node,
                                        const EList<SpliceSite>& spliceSites_db,
                                        uint32_t ref,
                                        uint32_t left,
                                        uint32_t right,
                                        EList<SpliceSite>& spliceSites) const
{
    assert(node != NULL);
    bool goleft = true, goright = true;
    if((node->key.ref() > ref) ||
       (node->key.ref() == ref && node->key.left() > right)) {
        goright = false;
    }
    
    if((node->key.ref() < ref) ||
       (node->key.ref() == ref && node->key.left() < left)) {
        goleft = false;
    }
    
    if(goleft && node->left != NULL) {
        getSpliceSites_recur(
                             node->left,
			     spliceSites_db,
                             ref,
                             left,
                             right,
                             spliceSites);
    }
    
    if(node->key.ref() == ref &&
       node->key.left() >= left && node->key.left() <= right) {
        assert_lt(node->payload, spliceSites_db.size());
#ifndef NDEBUG
        const SpliceSite& ss = spliceSites_db[node->payload];
        assert_eq(ss.ref(), node->key.ref());
        assert(ss.left() == node->key.left() ||
               ss.right() == node->key.left());
#endif
        spliceSites.push_back(spliceSites_db[node->payload]);
    }
    
    if(goright && node->right != NULL) {
        getSpliceSites_recur(
                             node->right,
			     spliceSites_db,
                             ref,
                             left,
                             right,
                             spliceSites);
    }
}

void SpliceSiteDB::print(ofstream& out)
{
    EList<SpliceSitePos> ss_list;
    for(uint64_t i = 0; i < _fwIndex.size(); i++) {
      for(uint64_t j = 0; j < _fwIndex[i].size(); j++) {
        assert(_fwIndex[i][j] != NULL);
        const Node *root = _fwIndex[i][j]->root();
        if(root != NULL) print_recur(root, out, ss_list);
      }
    }
    print_impl(out, ss_list);
}

void SpliceSiteDB::print_recur(
                               const RedBlackNode<SpliceSitePos, uint32_t> *node,
                               ofstream& out,
                               EList<SpliceSitePos>& ss_list)
{
    if(node == NULL) return;
    print_recur(node->left, out, ss_list);
    print_impl(out, ss_list, &(node->key));
    print_recur(node->right, out, ss_list);
}

void SpliceSiteDB::print_impl(
                              ofstream& out,
                              EList<SpliceSitePos>& ss_list,
                              const SpliceSitePos* ss)
{
    size_t i = 0;
    while(i < ss_list.size()) {
        const SpliceSitePos& tmp_ss = ss_list[i];
        bool do_print = true;
        if(ss != NULL) {
            if(tmp_ss.ref() == ss->ref()) {
                assert_leq(tmp_ss.left(), ss->left());
                if(ss->left() < tmp_ss.left() + 10) {
                    do_print = false;
                    if((int)ss->left() - (int)tmp_ss.left() == (int)ss->right() - (int)tmp_ss.right()) {
                        if(!tmp_ss.canonical() && ss->canonical()) {
                            ss_list.erase(i);
                            ss_list.push_back(*ss);
                        }
                        return;
                    }
                }
            }
        }
        
        if(!do_print) {
            i++;
            continue;
        }
        
        assert_lt(tmp_ss.ref(), _refnames.size());
        out << _refnames[tmp_ss.ref()] << "\t"
            << tmp_ss.left() << "\t"
            << tmp_ss.right() << "\t"
            << (tmp_ss.canonical() ? (tmp_ss.fw() ? "+" : "-") : ".") << endl;

        ss_list.erase(i);
    }
    
    if(ss != NULL) ss_list.push_back(*ss);
}

void SpliceSiteDB::read(ifstream& in, bool novel)
{
    _empty = false;
    assert_eq(_numRefs, _refnames.size());
    while(!in.eof()) {
        string refname;
        uint32_t left = 0, right = 0;
        char fw = 0;
        in >> refname >> left >> right >> fw;
        uint32_t ref = 0;
        for(; ref < _refnames.size(); ref++) {
            if(_refnames[ref] == refname) break;
        }
        if(ref >= _numRefs) continue;
        uint64_t ref_seg = left / ref_segment_size;
        assert_lt(ref, _spliceSites.size());
        assert_lt(ref_seg, _spliceSites[ref].size());
	
        _spliceSites[ref][ref_seg].expand();
        _spliceSites[ref][ref_seg].back().init(ref, left, right, fw == '+' || fw == '.', fw != '.');
        _spliceSites[ref][ref_seg].back()._fromfile = true;
        assert_gt(_spliceSites[ref][ref_seg].size(), 0);

        bool added = false;
        assert_lt(ref, _fwIndex.size());
        assert_lt(ref_seg, _fwIndex[ref].size());
        assert(_fwIndex[ref][ref_seg] != NULL);
        Node *cur = _fwIndex[ref][ref_seg]->add(pool(ref, ref_seg), _spliceSites[ref][ref_seg].back(), &added);
        assert(added);
        assert(cur != NULL);
        cur->payload = _spliceSites[ref][ref_seg].size() - 1;
        
        added = false;
        uint64_t ref_seg2 = right / ref_segment_size;
        SpliceSitePos rssp(ref, right, left, fw == '+' || fw == '.', fw != '.');
        assert_lt(ref, _bwIndex.size());
        assert_lt(ref_seg2, _bwIndex.size());
        assert(_bwIndex[ref][ref_seg2] != NULL);
        cur = _bwIndex[ref][ref_seg2]->add(pool(ref, ref_seg2), rssp, &added);
        assert(added);
        assert(cur != NULL);
        if(ref_seg != ref_seg2) {
            _spliceSites[ref][ref_seg2].expand();
            _spliceSites[ref][ref_seg2].back().init(ref, left, right, fw == '+' || fw == '.', fw != '.');
            _spliceSites[ref][ref_seg2].back()._fromfile = true;
        }
        cur->payload = _spliceSites[ref][ref_seg2].size() - 1;
    }
}

Pool& SpliceSiteDB::pool(uint64_t ref, uint64_t ref_seg) {
    assert_lt(ref, _numRefs);
    assert_lt(ref, _pool.size());
    assert_lt(ref_seg, _pool[ref].size());
    if(_pool[ref][ref_seg].size() <= 0 || _pool[ref][ref_seg].back()->full()) {
        _pool[ref][ref_seg].push_back(new Pool(1 << 18 /* 256KB */, 16 << 10 /* 16KB */, CA_CAT));
    }
    assert(_pool[ref][ref_seg].back() != NULL);
    return *_pool[ref][ref_seg].back();
}

#else

SpliceSiteDB::SpliceSiteDB(
                           const BitPairReference& refs,
                           const EList<string>& refnames,
                           bool threadSafe,
                           bool write,
                           bool read) :
_numRefs(refs.numRefs()),
_refnames(refnames),
_write(write),
_read(read),
_threadSafe(threadSafe),
_empty(true)
{
    assert_gt(_numRefs, 0);
    assert_eq(_numRefs, _refnames.size());
    for(uint64_t i = 0; i < _numRefs; i++) {
        _fwIndex.push_back(new RedBlack<SpliceSitePos, uint32_t>(16 << 10, CA_CAT));
        _bwIndex.push_back(new RedBlack<SpliceSitePos, uint32_t>(16 << 10, CA_CAT));
        _pool.expand();
        _spliceSites.expand();
        _mutex.push_back(MUTEX_T());
    }
    
    donorstr.resize(donor_exonic_len + donor_intronic_len);
    acceptorstr.resize(acceptor_intronic_len + acceptor_exonic_len);
}

SpliceSiteDB::~SpliceSiteDB() {
    assert_eq(_fwIndex.size(), _bwIndex.size());
    assert_eq(_fwIndex.size(), _pool.size());
    for(uint64_t i = 0; i < _numRefs; i++) {
        delete _fwIndex[i];
        delete _bwIndex[i];
        
        EList<Pool*>& pool = _pool[i];
        for(size_t j = 0; j < pool.size(); j++) {
            delete pool[j];
        }
    }
}

size_t SpliceSiteDB::size(uint64_t ref) const {
    if(!_read) return 0;
    
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    assert_lt(ref, _fwIndex.size());
    assert_eq(_fwIndex.size(), _bwIndex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    return _fwIndex.size();
}

bool SpliceSiteDB::empty(uint64_t ref) const {
    return size(ref) == 0;
}

bool SpliceSiteDB::addSpliceSite(
                                 const Read& rd,
                                 const AlnRes& rs,
                                 uint32_t minAnchorLen)
{
    if(!_write) return false;
    if(rs.trimmed5p(true) + rs.trimmed3p(true) > 0) return false;
    
    _empty = false;
    
    Coord coord = rs.refcoord();
    uint64_t ref = coord.ref();
    assert_lt(ref, _numRefs);
    const EList<Edit>& edits = rs.ned();
    if(!coord.orient()) {
        Edit::invertPoss(const_cast<EList<Edit>&>(edits), rd.length(), false);
    }
    SpliceSitePos ssp;
    uint32_t refoff = coord.off();
    uint32_t leftAnchorLen = 0, rightAnchorLen = 0;
    size_t eidx = 0;
    size_t last_eidx = 0;
    uint32_t mm = 0;
    for(size_t i = 0; i < rd.length(); i++, refoff++) {
        while(eidx < edits.size() && edits[eidx].pos == i) {
            if(edits[eidx].isReadGap()) {
                refoff++;
            } else if(edits[eidx].isRefGap()) {
                assert_gt(refoff, 0);
                refoff--;
            }
            if(edits[eidx].isGap() || edits[eidx].isMismatch()) mm++;
            if(edits[eidx].isSpliced()) {
                assert_gt(refoff, 0);
                if(ssp.inited()) {
                    assert(edits[last_eidx].isSpliced());
                    assert_lt(edits[last_eidx].pos, edits[eidx].pos);
                    rightAnchorLen = edits[eidx].pos - edits[last_eidx].pos;
                    uint32_t minLeftAnchorLen = minAnchorLen + mm * 2 + (edits[eidx].splDir == EDIT_SPL_UNKNOWN ? 6 : 0);
                    uint32_t mm2 = 0;
                    for(size_t j = eidx + 1; j < edits.size(); j++) {
                        if(edits[j].isGap() || edits[j].isMismatch()) mm2++;
                    }
                    uint32_t minRightAnchorLen = minAnchorLen + mm2 * 2 + (edits[eidx].splDir == EDIT_SPL_UNKNOWN ? 6 : 0);
                    if(leftAnchorLen >= minLeftAnchorLen && rightAnchorLen >= minRightAnchorLen) {
                        bool added = false;
                        assert_lt(ref, _mutex.size());
                        ThreadSafe t(&_mutex[ref], _threadSafe && _write);
                        assert_lt(ref, _fwIndex.size());
                        assert(_fwIndex[ref] != NULL);
                        Node *cur = _fwIndex[ref]->add(pool(ref), ssp, &added);
                        if(added) {
                            assert_lt(ref, _spliceSites.size());
                            _spliceSites[ref].expand();
                            _spliceSites[ref].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                            _spliceSites[ref].back()._readid = rd.rdid;
                            _spliceSites[ref].back()._leftext = leftAnchorLen;
                            _spliceSites[ref].back()._rightext = rightAnchorLen;
                            _spliceSites[ref].back()._numreads = 1;
                            assert(cur != NULL);
                            cur->payload = _spliceSites[ref].size() - 1;
                            
                            SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                            assert_lt(ref, _bwIndex.size());
                            assert(_bwIndex[ref] != NULL);
                            cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
                            assert(added);
                            assert(cur != NULL);
                            cur->payload = _spliceSites[ref].size() - 1;
                            assert_eq(_fwIndex[ref]->size(), _bwIndex[ref]->size());
                        } else {
                            assert(cur != NULL);
                            assert_lt(ref, _spliceSites.size());
                            assert_lt(cur->payload, _spliceSites[ref].size());
                            if(leftAnchorLen > _spliceSites[ref][cur->payload]._leftext) _spliceSites[ref][cur->payload]._leftext = leftAnchorLen;
                            if(rightAnchorLen > _spliceSites[ref][cur->payload]._rightext) _spliceSites[ref][cur->payload]._rightext = rightAnchorLen;
                            _spliceSites[ref][cur->payload]._numreads += 1;
                            if(rd.rdid < _spliceSites[ref][cur->payload]._readid) {
                                _spliceSites[ref][cur->payload]._readid = rd.rdid;
                            }
                        }
                    }
                    leftAnchorLen = rightAnchorLen;
                    rightAnchorLen = 0;
                } else {
                    leftAnchorLen = edits[eidx].pos;
                }
                bool fw = (edits[eidx].splDir == EDIT_SPL_FW || edits[eidx].splDir == EDIT_SPL_UNKNOWN);
                bool canonical = (edits[eidx].splDir != EDIT_SPL_UNKNOWN);
                ssp.init(coord.ref(), refoff - 1, refoff + edits[eidx].splLen, fw, canonical);
                refoff += edits[eidx].splLen;
                last_eidx = eidx;
            }
            eidx++;
        }
    }
    if(ssp.inited()) {
        assert(edits[last_eidx].isSpliced());
        assert_lt(edits[last_eidx].pos, rd.length());
        rightAnchorLen = rd.length() - edits[last_eidx].pos;
        uint32_t minLeftAnchorLen = minAnchorLen + mm * 2 + (edits[last_eidx].splDir == EDIT_SPL_UNKNOWN ? 6 : 0);
        uint32_t mm2 = 0;
        for(size_t j = last_eidx + 1; j < edits.size(); j++) {
            if(edits[j].isGap() || edits[j].isMismatch()) mm2++;
        }
        uint32_t minRightAnchorLen = minAnchorLen + mm2 * 2 + (edits[last_eidx].splDir == EDIT_SPL_UNKNOWN ? 6 : 0);
        if(leftAnchorLen >= minLeftAnchorLen && rightAnchorLen >= minRightAnchorLen) {
            bool added = false;
            assert_lt(ref, _mutex.size());
            ThreadSafe t(&_mutex[ref], _threadSafe && _write);
            assert_lt(ref, _fwIndex.size());
            assert(_fwIndex[ref] != NULL);
            Node *cur = _fwIndex[ref]->add(pool(ref), ssp, &added);
            if(added) {
                assert_lt(ref, _spliceSites.size());
                _spliceSites[ref].expand();
                _spliceSites[ref].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.fw(), ssp.canonical());
                _spliceSites[ref].back()._readid = rd.rdid;
                _spliceSites[ref].back()._leftext = leftAnchorLen;
                _spliceSites[ref].back()._rightext = rightAnchorLen;
                _spliceSites[ref].back()._numreads = 1;
                assert(cur != NULL);
                cur->payload = _spliceSites[ref].size() - 1;
                
                SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.fw(), ssp.canonical());
                assert_lt(ref, _bwIndex.size());
                assert(_bwIndex[ref] != NULL);
                cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
                assert(added);
                assert(cur != NULL);
                cur->payload = _spliceSites[ref].size() - 1;
                assert_eq(_fwIndex[ref]->size(), _bwIndex[ref]->size());
            } else {
                assert(cur != NULL);
                assert_lt(ref, _spliceSites.size());
                assert_lt(cur->payload, _spliceSites[ref].size());
                if(leftAnchorLen > _spliceSites[ref][cur->payload]._leftext) _spliceSites[ref][cur->payload]._leftext = leftAnchorLen;
                if(rightAnchorLen > _spliceSites[ref][cur->payload]._rightext) _spliceSites[ref][cur->payload]._rightext = rightAnchorLen;
                _spliceSites[ref][cur->payload]._numreads += 1;
                if(rd.rdid < _spliceSites[ref][cur->payload]._readid) {
                    _spliceSites[ref][cur->payload]._readid = rd.rdid;
                }
            }
        }
    }
    if(!coord.orient()) {
        Edit::invertPoss(const_cast<EList<Edit>&>(edits), rd.length(), false);
    }
    
    return true;
}

bool SpliceSiteDB::getSpliceSite(SpliceSite& ss) const
{
    if(!_read) return false;
    
    uint64_t ref = ss.ref();
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    
    assert_lt(ref, _fwIndex.size());
    assert(_fwIndex[ref] != NULL);
    const Node *cur = _fwIndex[ref]->lookup(ss);
    if(cur == NULL) return false;
    assert(cur != NULL);
    assert_lt(ref, _spliceSites.size());
    ss = _spliceSites[ref][cur->payload];
    return true;
}

void SpliceSiteDB::getLeftSpliceSites(uint32_t ref, uint32_t left, uint32_t range, EList<SpliceSite>& spliceSites) const
{
    if(!_read) return;
    
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    assert_gt(range, 0);
    assert_geq(left + 1, range);
    assert_lt(ref, _bwIndex.size());
    assert(_bwIndex[ref] != NULL);
    const Node *cur = _bwIndex[ref]->root();
    if(cur != NULL) getSpliceSites_recur(cur, ref, left + 1 - range, left, spliceSites);
}

void SpliceSiteDB::getRightSpliceSites(uint32_t ref, uint32_t right, uint32_t range, EList<SpliceSite>& spliceSites) const
{
    if(!_read) return;
    
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    assert_gt(range, 0);
    assert_gt(right + range, range);
    assert_lt(ref, _fwIndex.size());
    assert(_fwIndex[ref] != NULL);
    const Node *cur = _fwIndex[ref]->root();
    if(cur != NULL) getSpliceSites_recur(cur, ref, right, right + range - 1, spliceSites);
    
}

void SpliceSiteDB::getSpliceSites_recur(
                                        const RedBlackNode<SpliceSitePos, uint32_t> *node,
                                        uint32_t ref,
                                        uint32_t left,
                                        uint32_t right,
                                        EList<SpliceSite>& spliceSites) const
{
    assert(node != NULL);
    bool goleft = true, goright = true;
    if((node->key.ref() > ref) ||
       (node->key.ref() == ref && node->key.left() > right)) {
        goright = false;
    }
    
    if((node->key.ref() < ref) ||
       (node->key.ref() == ref && node->key.left() < left)) {
        goleft = false;
    }
    
    if(goleft && node->left != NULL) {
        getSpliceSites_recur(
                             node->left,
                             ref,
                             left,
                             right,
                             spliceSites);
    }
    
    if(node->key.ref() == ref &&
       node->key.left() >= left && node->key.left() <= right) {
        assert_lt(ref, _spliceSites.size());
        assert_lt(node->payload, _spliceSites[ref].size());
#ifndef NDEBUG
        const SpliceSite& ss = _spliceSites[ref][node->payload];
        assert_eq(ss.ref(), node->key.ref());
        assert(ss.left() == node->key.left() ||
               ss.right() == node->key.left());
#endif
        spliceSites.push_back(_spliceSites[ref][node->payload]);
    }
    
    if(goright && node->right != NULL) {
        getSpliceSites_recur(
                             node->right,
                             ref,
                             left,
                             right,
                             spliceSites);
    }
}

void calculate_splicesite_read_dist_impl(const RedBlackNode<SpliceSitePos, uint32_t> *node,
                                         const EList<SpliceSite> &spliceSites,
                                         EList<int64_t>& splicesite_read_dist) {
    if(node == NULL) return;
    calculate_splicesite_read_dist_impl(node->left, spliceSites, splicesite_read_dist);
    assert_lt(node->payload, spliceSites.size());
    const SpliceSite& ss = spliceSites[node->payload];
    if(ss.numreads() < splicesite_read_dist.size())
        splicesite_read_dist[ss.numreads()] += 1;
    else
        splicesite_read_dist.back() += 1;
    calculate_splicesite_read_dist_impl(node->right, spliceSites, splicesite_read_dist);
}

uint32_t calculate_splicesite_read_dist(const EList<RedBlack<SpliceSitePos, uint32_t>* >& fwIndex,
                                    const ELList<SpliceSite> &spliceSites,
                                    EList<int64_t>& splicesite_read_dist) {
    for(size_t i = 0; i < fwIndex.size(); i++) {
        assert(fwIndex[i] != NULL);
        const RedBlackNode<SpliceSitePos, uint32_t> *root = fwIndex[i]->root();
        assert_lt(i, spliceSites.size());
        if(root != NULL) calculate_splicesite_read_dist_impl(root, spliceSites[i], splicesite_read_dist);
    }
    
    for(size_t i = 1; i < splicesite_read_dist.size(); i++) {
        splicesite_read_dist[i] += splicesite_read_dist[i-1];
    }
    
    for(size_t i = 0; i < splicesite_read_dist.size(); i++) {
        float cmf_i = float(splicesite_read_dist[i]) / splicesite_read_dist.back();
        if(cmf_i > 0.7)
            return i;
    }
    
    return 0;
}

void SpliceSiteDB::print(ofstream& out)
{
    EList<int64_t> splicesite_read_dist;
    for(size_t i = 0; i < 100; i++) {
        splicesite_read_dist.push_back(0);
    }
    uint32_t numreads_cutoff = calculate_splicesite_read_dist(_fwIndex, _spliceSites, splicesite_read_dist);
    
    EList<SpliceSite> ss_list;
    for(size_t i = 0; i < _fwIndex.size(); i++) {
        assert(_fwIndex[i] != NULL);
        const Node *root = _fwIndex[i]->root();
        if(root != NULL) print_recur(root, out, numreads_cutoff, ss_list);
    }
    print_impl(out, ss_list);
}

void SpliceSiteDB::print_recur(
                               const RedBlackNode<SpliceSitePos, uint32_t> *node,
                               ofstream& out,
                               const uint32_t numreads_cutoff,
                               EList<SpliceSite>& ss_list)
{
    if(node == NULL) return;
    print_recur(node->left, out, numreads_cutoff, ss_list);
    const SpliceSitePos& ssp = node->key;
    assert_lt(ssp.ref(), _spliceSites.size());
    assert_lt(node->payload, _spliceSites[ssp.ref()].size());
    const SpliceSite& ss = _spliceSites[ssp.ref()][node->payload];
    if(ss.numreads() >= numreads_cutoff) print_impl(out, ss_list, &ss);
    print_recur(node->right, out, numreads_cutoff, ss_list);
}

void SpliceSiteDB::print_impl(
                              ofstream& out,
                              EList<SpliceSite>& ss_list,
                              const SpliceSite* ss)
{
    size_t i = 0;
    while(i < ss_list.size()) {
        const SpliceSite& tmp_ss = ss_list[i];
        bool do_print = true;
        if(ss != NULL) {
            if(tmp_ss.ref() == ss->ref()) {
                assert_leq(tmp_ss.left(), ss->left());
                if(ss->left() < tmp_ss.left() + 10) {
                    do_print = false;
                    if(abs(((int)ss->left() - (int)tmp_ss.left()) - ((int)ss->right() - (int)tmp_ss.right())) <= 10) {
                        if(tmp_ss.numreads() < ss->numreads()) {
                            ss_list.erase(i);
                            ss_list.push_back(*ss);
                        }
                        return;
                    }
                }
            }
        }
        
        if(!do_print) {
            i++;
            continue;
        }
        
        assert_lt(tmp_ss.ref(), _refnames.size());
        out << _refnames[tmp_ss.ref()] << "\t"
        << tmp_ss.left() << "\t"
        << tmp_ss.right() << "\t"
        << (tmp_ss.canonical() ? (tmp_ss.fw() ? "+" : "-") : ".") << endl;
        
        ss_list.erase(i);
    }
    
    if(ss != NULL) ss_list.push_back(*ss);
}

void SpliceSiteDB::read(ifstream& in, bool novel)
{
    _empty = false;
    assert_eq(_numRefs, _refnames.size());
    while(!in.eof()) {
        string refname;
        uint32_t left = 0, right = 0;
        char fw = 0;
        in >> refname >> left >> right >> fw;
        uint32_t ref = 0;
        for(; ref < _refnames.size(); ref++) {
            if(_refnames[ref] == refname) break;
        }
        if(ref >= _numRefs) continue;
        assert_lt(ref, _spliceSites.size());
        _spliceSites[ref].expand();
        _spliceSites[ref].back().init(ref, left, right, fw == '+' || fw == '.', fw != '.');
        _spliceSites[ref].back()._fromfile = true;
        assert_gt(_spliceSites[ref].size(), 0);
        
        bool added = false;
        assert_lt(ref, _fwIndex.size());
        assert(_fwIndex[ref] != NULL);
        Node *cur = _fwIndex[ref]->add(pool(ref), _spliceSites[ref].back(), &added);
        assert(added);
        assert(cur != NULL);
        cur->payload = _spliceSites[ref].size() - 1;
        
        added = false;
        SpliceSitePos rssp(ref, right, left, fw == '+' || fw == '.', fw != '.');
        assert_lt(ref, _bwIndex.size());
        assert(_bwIndex[ref] != NULL);
        cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
        assert(added);
        assert(cur != NULL);
        cur->payload = _spliceSites[ref].size() - 1;
    }
}

Pool& SpliceSiteDB::pool(uint64_t ref) {
    assert_lt(ref, _numRefs);
    assert_lt(ref, _pool.size());
    EList<Pool*>& pool = _pool[ref];
    if(pool.size() <= 0 || pool.back()->full()) {
        pool.push_back(new Pool(1 << 20 /* 1MB */, 16 << 10 /* 16KB */, CA_CAT));
    }
    assert(pool.back() != NULL);
    return *pool.back();
}

#endif

float SpliceSiteDB::probscore(
                              int64_t donor_seq,
                              int64_t acceptor_seq)
{
    float probscore = 0.0f;
#if defined(NEW_PROB_MODEL)
    float donor_probscore = 0.0f;
    assert_leq(donor_seq, 0x3ffff);
    int64_t donor_exonic_seq = (donor_seq >> 4) & (~0xff);
    int64_t donor_intronic_seq = donor_seq & 0xff;
    int64_t donor_rest_seq = donor_exonic_seq | donor_intronic_seq;
    int donor_seq3 = (donor_seq >> 10) & 0x3;
    int donor_seq4 = (donor_seq >> 8) & 0x3;
    donor_probscore = donor_cons1[donor_seq3] * donor_cons2[donor_seq4] / (background_bp_prob[donor_seq3] * background_bp_prob[donor_seq4]) * donor_me2x5[donor_rest_seq];
    
    float acceptor_probscore = 0.0f;
    assert_leq(acceptor_seq, 0x3fffffffffff);
    int64_t acceptor_intronic_seq = (acceptor_seq >> 4) & (~0x3f);
    int64_t acceptor_exonic_seq = acceptor_seq & 0x3f;
    int64_t acceptor_rest_seq = acceptor_intronic_seq | acceptor_exonic_seq;
    int acceptor_seq18 = (acceptor_seq >> 8) & 0x3;
    int acceptor_seq19 = (acceptor_seq >> 6) & 0x3;
    acceptor_probscore = acceptor_cons1[acceptor_seq18] * acceptor_cons2[acceptor_seq19] / (background_bp_prob[acceptor_seq18] * background_bp_prob[acceptor_seq19]);
    
    int64_t acceptor_seq1 = acceptor_rest_seq >> 28 & 0x3fff; // [0, 7]
    acceptor_probscore *= acceptor_me2x3acc1[acceptor_seq1];
    int64_t acceptor_seq2 = (acceptor_rest_seq >> 14) & 0x3fff; // [7, 7]
    acceptor_probscore *= acceptor_me2x3acc2[acceptor_seq2];
    int64_t acceptor_seq3 = acceptor_rest_seq & 0x3fff; // [14, 7]
    acceptor_probscore *= acceptor_me2x3acc3[acceptor_seq3];
    int64_t acceptor_seq4 = (acceptor_rest_seq >> 20) & 0x3fff; // [4, 7]
    acceptor_probscore *= acceptor_me2x3acc4[acceptor_seq4];
    int64_t acceptor_seq5 = (acceptor_rest_seq >> 6) & 0x3fff; // [11, 7]
    acceptor_probscore *= acceptor_me2x3acc5[acceptor_seq5];
    int64_t acceptor_seq6 = acceptor_seq1 & 0x3f; // [4, 3]
    acceptor_probscore /= acceptor_me2x3acc6[acceptor_seq6];
    int64_t acceptor_seq7 = acceptor_seq4 & 0xff; // [7, 4]
    acceptor_probscore /= acceptor_me2x3acc7[acceptor_seq7];
    int64_t acceptor_seq8 = acceptor_seq2 & 0x3f; // [11, 3]
    acceptor_probscore /= acceptor_me2x3acc8[acceptor_seq8];
    int64_t acceptor_seq9 = acceptor_seq5 & 0xff; // [14, 4]
    acceptor_probscore /= acceptor_me2x3acc9[acceptor_seq9];

    donor_probscore /= (1.0f + donor_probscore);
    acceptor_probscore /= (1.0f + acceptor_probscore);
    probscore = (donor_probscore + acceptor_probscore) / 2.0;
    
#else
    assert_lt(donor_seq, (int)(1 << (donor_len << 1)));
    probscore = donor_prob_sum[donor_seq];
    
    int acceptor_seq1 = acceptor_seq >> (acceptor_len2 << 1);
    assert_lt(acceptor_seq1, (int)(1 << (acceptor_len1 << 1)));
    probscore *= acceptor_prob_sum1[acceptor_seq1];
    
    int acceptor_seq2 = acceptor_seq % (1 << (acceptor_len2 << 1));
    probscore *= acceptor_prob_sum2[acceptor_seq2];
    
    probscore = 1.0 / (1.0 + probscore);
#endif
    return probscore;
}

