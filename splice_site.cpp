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
        << s.right() << "\t";
    if(s.splDir() == SPL_FW || s.splDir() == SPL_SEMI_FW) {
        out << "+";
    } else if(s.splDir() == SPL_RC || s.splDir() == SPL_SEMI_RC) {
        out << "-";
    } else {
        out << ".";
    }
    out << endl;
	return out;
}

SpliceSiteDB::SpliceSiteDB(
                           const BitPairReference& refs,
                           const EList<string>& refnames,
                           bool threadSafe,
                           bool write,
                           bool read) :
_numRefs(refs.numRefs()),
_write(write),
_read(read),
_threadSafe(threadSafe),
_empty(true)
{
    for(size_t r = 0; r < refnames.size(); r++) {
        const string& refname = refnames[r];
        _refnames.expand();
        size_t i = 0;
        for(; i < refname.size(); i++) {
            if(isspace(refname[i])) {
                break;
            }
        }
         _refnames.back() = refname.substr(0, i);
    }
    
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
    // daehwan - for debugging purposes
    uint32_t editdist = 0;
    for(size_t i = 0; i < edits.size(); i++) {
        const Edit& edit = edits[i];
        if(edit.isGap() || edit.isMismatch()) editdist++;
    }
        
    SpliceSitePos ssp;
    uint32_t refoff = (uint32_t)coord.off();
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
                    uint32_t minLeftAnchorLen = minAnchorLen + mm * 2 + (edits[eidx].splDir == SPL_UNKNOWN ? 6 : 0);
                    uint32_t mm2 = 0;
                    for(size_t j = eidx + 1; j < edits.size(); j++) {
                        if(edits[j].isGap() || edits[j].isMismatch()) mm2++;
                    }
                    uint32_t minRightAnchorLen = minAnchorLen + mm2 * 2 + (edits[eidx].splDir == SPL_UNKNOWN ? 6 : 0);
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
                            _spliceSites[ref].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.splDir());
                            _spliceSites[ref].back()._readid = rd.rdid;
                            _spliceSites[ref].back()._leftext = leftAnchorLen;
                            _spliceSites[ref].back()._rightext = rightAnchorLen;
                            _spliceSites[ref].back()._editdist = editdist;
                            _spliceSites[ref].back()._numreads = 1;
                            assert(cur != NULL);
                            cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
                            
                            SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.splDir());
                            assert_lt(ref, _bwIndex.size());
                            assert(_bwIndex[ref] != NULL);
                            cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
                            assert(added);
                            assert(cur != NULL);
                            cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
                            assert_eq(_fwIndex[ref]->size(), _bwIndex[ref]->size());
                        } else {
                            assert(cur != NULL);
                            assert_lt(ref, _spliceSites.size());
                            assert_lt(cur->payload, _spliceSites[ref].size());
                            if(leftAnchorLen > _spliceSites[ref][cur->payload]._leftext) _spliceSites[ref][cur->payload]._leftext = leftAnchorLen;
                            if(rightAnchorLen > _spliceSites[ref][cur->payload]._rightext) _spliceSites[ref][cur->payload]._rightext = rightAnchorLen;
                            if(editdist < _spliceSites[ref][cur->payload]._editdist) _spliceSites[ref][cur->payload]._editdist = editdist;
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
                ssp.init((uint32_t)coord.ref(), refoff - 1, refoff + edits[eidx].splLen, edits[eidx].splDir);
                refoff += edits[eidx].splLen;
                last_eidx = eidx;
            }
            eidx++;
        }
    }
    if(ssp.inited()) {
        assert(edits[last_eidx].isSpliced());
        assert_lt(edits[last_eidx].pos, rd.length());
        rightAnchorLen = (uint32_t)(rd.length() - edits[last_eidx].pos);
        uint32_t minLeftAnchorLen = minAnchorLen + mm * 2 + (edits[last_eidx].splDir == SPL_UNKNOWN ? 6 : 0);
        uint32_t mm2 = 0;
        for(size_t j = last_eidx + 1; j < edits.size(); j++) {
            if(edits[j].isGap() || edits[j].isMismatch()) mm2++;
        }
        uint32_t minRightAnchorLen = minAnchorLen + mm2 * 2 + (edits[last_eidx].splDir == SPL_UNKNOWN ? 6 : 0);
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
                _spliceSites[ref].back().init(ssp.ref(), ssp.left(), ssp.right(), ssp.splDir());
                _spliceSites[ref].back()._readid = rd.rdid;
                _spliceSites[ref].back()._leftext = leftAnchorLen;
                _spliceSites[ref].back()._rightext = rightAnchorLen;
                _spliceSites[ref].back()._editdist = editdist;
                _spliceSites[ref].back()._numreads = 1;
                assert(cur != NULL);
                cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
                
                SpliceSitePos rssp(ssp.ref(), ssp.right(), ssp.left(), ssp.splDir());
                assert_lt(ref, _bwIndex.size());
                assert(_bwIndex[ref] != NULL);
                cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
                assert(added);
                assert(cur != NULL);
                cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
                assert_eq(_fwIndex[ref]->size(), _bwIndex[ref]->size());
            } else {
                assert(cur != NULL);
                assert_lt(ref, _spliceSites.size());
                assert_lt(cur->payload, _spliceSites[ref].size());
                if(leftAnchorLen > _spliceSites[ref][cur->payload]._leftext) _spliceSites[ref][cur->payload]._leftext = leftAnchorLen;
                if(rightAnchorLen > _spliceSites[ref][cur->payload]._rightext) _spliceSites[ref][cur->payload]._rightext = rightAnchorLen;
                if(editdist < _spliceSites[ref][cur->payload]._editdist) _spliceSites[ref][cur->payload]._editdist = editdist;
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
    if(cur != NULL) getSpliceSites_recur(cur, left + 1 - range, left, spliceSites);
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
    if(cur != NULL) getSpliceSites_recur(cur, right, right + range - 1, spliceSites);
    
}

void SpliceSiteDB::getSpliceSites_recur(
                                        const RedBlackNode<SpliceSitePos, uint32_t> *node,
                                        uint32_t left,
                                        uint32_t right,
                                        EList<SpliceSite>& spliceSites) const
{
    assert(node != NULL);
    if(node->key.left() >= left && node->left != NULL) {
        getSpliceSites_recur(
                             node->left,
                             left,
                             right,
                             spliceSites);
    }
    
    if(node->key.left() >= left && node->key.left() <= right) {
        uint32_t ref = node->key.ref();
        assert_lt(ref, _spliceSites.size());
        assert_lt(node->payload, _spliceSites[ref].size());
        ASSERT_ONLY(const SpliceSite& ss = _spliceSites[ref][node->payload]);
        assert_eq(ss.ref(), node->key.ref());
        assert(ss.left() == node->key.left() ||
               ss.right() == node->key.left());
         spliceSites.push_back(_spliceSites[ref][node->payload]);
    }
    
    if(node->key.left() <= right && node->right != NULL) {
        getSpliceSites_recur(
                             node->right,
                             left,
                             right,
                             spliceSites);
    }
}

bool SpliceSiteDB::hasSpliceSites(
                                  uint32_t ref,
                                  uint32_t left1,
                                  uint32_t right1,
                                  uint32_t left2,
                                  uint32_t right2,
                                  bool includeNovel) const
{
    if(!_read) return false;
    
    assert_lt(ref, _numRefs);
    assert_lt(ref, _mutex.size());
    ThreadSafe t(const_cast<MUTEX_T*>(&_mutex[ref]), _threadSafe && _write);
    
    if(left1 < right1) {
        assert_lt(ref, _bwIndex.size());
        assert(_bwIndex[ref] != NULL);
        const Node *cur = _bwIndex[ref]->root();
        if(cur != NULL) {
            if(hasSpliceSites_recur(cur, left1, right1, includeNovel))
                return true;
        }
    }

    if(left2 < right2) {
        assert_lt(ref, _fwIndex.size());
        assert(_fwIndex[ref] != NULL);
        const Node *cur = _fwIndex[ref]->root();
        if(cur != NULL) {
            return hasSpliceSites_recur(cur, left2, right2, includeNovel);
        }
    }
    return false;
}

bool SpliceSiteDB::hasSpliceSites_recur(
                                        const RedBlackNode<SpliceSitePos, uint32_t> *node,
                                        uint32_t left,
                                        uint32_t right,
                                        bool includeNovel) const
{
    assert(node != NULL);
    if(node->key.left() >= left && node->key.left() <= right) {
        uint32_t ref = node->key.ref();
        assert_lt(ref, _spliceSites.size());
        assert_lt(node->payload, _spliceSites[ref].size());
        const SpliceSite& ss = _spliceSites[ref][node->payload];
        if(includeNovel || ss._known)
            return true;
    }
    
    if(node->key.left() >= left && node->left != NULL) {
        if(hasSpliceSites_recur(
                                node->left,
                                left,
                                right,
                                includeNovel))
            return true;
    }
    
    if(node->key.left() <= right && node->right != NULL) {
        if(hasSpliceSites_recur(
                                node->right,
                                left,
                                right,
                                includeNovel))
            return true;
    }
    
    return false;
}

bool SpliceSiteDB::insideExon(
                              uint32_t ref,
                              uint32_t left,
                              uint32_t right) const
{
    if(_exons.empty()) return false;
    assert_lt(ref, _numRefs);
    assert_lt(left, right);
    
    Exon e(ref, left + 1, 0, true);
    size_t i = _exons.bsearchLoBound(e);
    for(; i > 0; i--) {
        const Exon& e = _exons[i-1];
        if(e.right() < left) break;
        if(e.left() <= left && right <= e.right())
            return true;
    }
    return false;
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
            return (uint32_t)i;
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
    size_t numsplicesites = 0;
    for(size_t i = 0; i < _spliceSites.size(); i++) {
        numsplicesites += _spliceSites[i].size();
    }
    uint32_t numreads_cutoff2 = (uint32_t)(numsplicesites / 100000);
    
    EList<SpliceSite> ss_list;
    for(size_t i = 0; i < _fwIndex.size(); i++) {
        assert(_fwIndex[i] != NULL);
        const Node *root = _fwIndex[i]->root();
        if(root != NULL) print_recur(root, out, numreads_cutoff, numreads_cutoff2, ss_list);
    }
    print_impl(out, ss_list);
}

void SpliceSiteDB::print_recur(
                               const RedBlackNode<SpliceSitePos, uint32_t> *node,
                               ofstream& out,
                               const uint32_t numreads_cutoff,
                               const uint32_t numreads_cutoff2,
                               EList<SpliceSite>& ss_list)
{
    if(node == NULL) return;
    print_recur(node->left, out, numreads_cutoff, numreads_cutoff2, ss_list);
    const SpliceSitePos& ssp = node->key;
    assert_lt(ssp.ref(), _spliceSites.size());
    assert_lt(node->payload, _spliceSites[ssp.ref()].size());
    const SpliceSite& ss = _spliceSites[ssp.ref()][node->payload];
    if(ss.numreads() >= numreads_cutoff ||
       (ss.editdist() == 0 && ss.numreads() >= numreads_cutoff2)) print_impl(out, ss_list, &ss);
    print_recur(node->right, out, numreads_cutoff, numreads_cutoff2, ss_list);
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
        << tmp_ss.right() << "\t";
        if(tmp_ss.splDir() == SPL_FW || tmp_ss.splDir() == SPL_SEMI_FW) {
            out << "+";
        } else if(tmp_ss.splDir() == SPL_RC || tmp_ss.splDir() == SPL_SEMI_RC) {
            out << "-";
        } else {
            out << ".";
        }
        out << endl;
        ss_list.erase(i);
    }
    
    if(ss != NULL) ss_list.push_back(*ss);
}

void SpliceSiteDB::read(const GFM<TIndexOffU>& gfm, const EList<ALT<TIndexOffU> >& alts)
{
    EList<Exon> exons;
    _empty = false;
    assert_eq(_numRefs, _refnames.size());
    for(size_t i = 0; i < alts.size(); i++) {
        const ALT<TIndexOffU>& alt = alts[i];
        if(!alt.splicesite() && !alt.exon()) continue;
        if(alt.left > alt.right) continue;
        TIndexOffU ref = 0, left = 0, tlen = 0;
        char fw = alt.fw;
        bool straddled2 = false;
        gfm.joinedToTextOff(
                            1,
                            alt.left,
                            ref,
                            left,
                            tlen,
                            true,         // reject straddlers?
                            straddled2);  // straddled?
        assert_lt(ref, _spliceSites.size());
        TIndexOffU right = left + (alt.right - alt.left);
        if(alt.splicesite()) {
            left -= 1; right += 1;
            _spliceSites[ref].expand();
            _spliceSites[ref].back().init(ref,
                                          left,
                                          right,
                                          fw ? SPL_FW : SPL_RC,
                                          alt.exon(),
                                          true,   // from file?
                                          true);  // known splice site?
            assert_gt(_spliceSites[ref].size(), 0);
            bool added = false;
            assert_lt(ref, _fwIndex.size());
            assert(_fwIndex[ref] != NULL);
            Node *cur = _fwIndex[ref]->add(pool(ref), _spliceSites[ref].back(), &added);
            if(!added) {
                _spliceSites[ref].pop_back();
                continue;
            }
            assert(added);
            assert(cur != NULL);
            cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
            
            added = false;
            SpliceSitePos rssp(ref,
                               right,
                               left,
                               fw ? SPL_FW : SPL_RC);
            assert_lt(ref, _bwIndex.size());
            assert(_bwIndex[ref] != NULL);
            cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
            assert(added);
            assert(cur != NULL);
            cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
        } else {
            assert(alt.exon());
            // Given some relaxation
            if(left >= 10) left -= 10;
            else           left = 0;
            if(right + 10 < tlen) right += 10;
            else                  right = tlen - 1;
            exons.expand();
            exons.back().init(ref, left, right, fw == '+' ? SPL_FW : SPL_RC);
        }
    }
    if(exons.size() > 0) {
        _exons.resizeExact(exons.size()); _exons.clear();
        _exons.push_back_array(exons.begin(), exons.size());
        _exons.sort();
    }
}

void SpliceSiteDB::read(ifstream& in, bool known)
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
        _spliceSites[ref].back().init(ref,
                                      left,
                                      right,
                                      fw == '+' ? SPL_FW : SPL_RC,
                                      false,  // exon?
                                      true,   // from file?
                                      known); // known splice site?
        assert_gt(_spliceSites[ref].size(), 0);
        
        bool added = false;
        assert_lt(ref, _fwIndex.size());
        assert(_fwIndex[ref] != NULL);
        Node *cur = _fwIndex[ref]->add(pool(ref), _spliceSites[ref].back(), &added);
        if(!added) {
            _spliceSites[ref].pop_back();
            continue;
        }
        
        assert(cur != NULL);
        cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
        added = false;
        SpliceSitePos rssp(ref,
                           right,
                           left,
                           fw == '+' ? SPL_FW : SPL_RC);
        assert_lt(ref, _bwIndex.size());
        assert(_bwIndex[ref] != NULL);
        cur = _bwIndex[ref]->add(pool(ref), rssp, &added);
        assert(added);
        assert(cur != NULL);
        cur->payload = (uint32_t)_spliceSites[ref].size() - 1;
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
    
    int acceptor_seq1 = (int)(acceptor_seq >> (acceptor_len2 << 1));
    assert_lt(acceptor_seq1, (int)(1 << (acceptor_len1 << 1)));
    probscore *= acceptor_prob_sum1[acceptor_seq1];
    
    int acceptor_seq2 = acceptor_seq % (1 << (acceptor_len2 << 1));
    probscore *= acceptor_prob_sum2[acceptor_seq2];
    
    probscore = 1.0 / (1.0 + probscore);
#endif
    return probscore;
}

