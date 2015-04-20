#ifndef _RLCSA_DOCARRAY_H
#define _RLCSA_DOCARRAY_H

#include <stack>
#include <vector>

#include "rlcsa.h"


namespace CSA
{

struct STNode
{
  uint      string_depth;
  pair_type range;
  STNode*   parent;
  STNode*   child;
  STNode*   sibling;
  STNode*   next;

  usint     stored_documents;
  uint      id;
  bool      contains_all;

  std::vector<uint>* docs;

  STNode(uint lcp, pair_type sa_range);
  ~STNode();

  void addChild(STNode* node);
  void addSibling(STNode* node);

  void deleteChildren();

  void addLeaves();   // Adds leaves to the sparse suffix tree when required.
  STNode* addLeaf(STNode* left, STNode* right, usint pos);

  bool verifyTree();  // Call this to ensure that the tree is correct.

  // Removes this node from the tree, replacing it by its children.
  void remove();

  void determineSize(uint& nodes, uint& leaves);
  void computeStoredDocuments(usint* documents);
  void containsAllDocuments();  // Use this to tell that the node contains all documents.
  void storeThisSet();  // Use this to tell that the set contained in the node will be stored.
  void setNext();
};

std::ostream& operator<<(std::ostream& stream, const STNode& data);

class DocArray
{
  public:
    const static usint PARALLEL_SORT_THRESHOLD = 4096;

    explicit DocArray(const RLCSA& _rlcsa); // This DocArray can only be used for direct listing. Nothing else works.
    DocArray(STNode* root, const RLCSA& _rlcsa);
    DocArray(const RLCSA& _rlcsa, const std::string& base_name, bool load_grammar = true);
    ~DocArray();

    void readRules(const std::string& name_prefix, bool print = false);

    void writeTo(const std::string& base_name) const;
    usint reportSize(bool print = false) const;

    std::vector<usint>* listDocuments(const std::string& pattern) const;
    std::vector<usint>* listDocuments(pair_type range) const;

    usint count(const std::string& pattern) const;
    usint count(pair_type sa_range) const;

    // These versions run-length encode the answers.
    std::vector<pair_type>* listDocumentsRLE(const std::string& pattern) const;
    std::vector<pair_type>* listDocumentsRLE(pair_type range) const;

    // This version uses RLCSA directly.
    std::vector<usint>* directListing(pair_type range) const;
    std::vector<pair_type>* directListingRLE(pair_type range) const;

    inline bool isOk() const { return this->ok; }
    inline bool hasGrammar() const { return this->has_grammar; }
    inline bool usesRLE() const { return this->uses_rle; }
    inline usint getSize() const { return this->rlcsa.getSize(); }

    inline usint getNumberOfNodes() const { return this->getNumberOfLeaves() + this->getNumberOfInternalNodes(); }
    inline usint getNumberOfLeaves() const { return this->leaf_ranges->getNumberOfItems(); }
    inline usint getNumberOfInternalNodes() const { return this->first_children->getNumberOfItems(); }

    inline usint maxInteger() const { return this->getNumberOfDocuments() + this->getNumberOfRules(); }
    inline usint getNumberOfDocuments() const { return this->rlcsa.getNumberOfSequences(); }
    inline usint getNumberOfRules() const { return this->rule_borders->getNumberOfItems(); }

//--------------------------------------------------------------------------

    /*
      In document graph, nodes 1 to ndoc are document ids and ndoc + 1 to ndoc + |nodes|
      are node/block ids.
    */

    inline bool nodeIsDoc(usint node_id) const
    {
      return (node_id > 0 && node_id <= this->getNumberOfDocuments());
    }

    inline bool nodeIsBlock(usint node_id) const
    {
      return (node_id > this->getNumberOfDocuments() &&
              node_id <= this->getNumberOfDocuments() + this->getNumberOfNodes());
    }

    inline usint nodeToDoc(usint node_id) const { return node_id - 1; }
    inline usint nodeToBlock(usint node_id) const { return node_id - this->getNumberOfDocuments() - 1; }

//--------------------------------------------------------------------------

  private:
    const static usint RANGE_BLOCK_SIZE = 32;
    const static usint PARENT_BLOCK_SIZE = 32;
    const static usint RULE_BLOCK_SIZE = 32;
    const static usint BLOCK_BLOCK_SIZE = 32;

    const RLCSA&    rlcsa;

    DeltaVector*    leaf_ranges;
    SuccinctVector* first_children;
    ReadBuffer*     parents;
    ReadBuffer*     next_leaves;  // this->getNumberOfLeaves() denotes that there is no next leaf.

    SuccinctVector* rule_borders;
    ReadBuffer*     rules;        // this->getNumberOfDocs() means all documents.

    SuccinctVector* block_borders;
    ReadBuffer*     blocks;       // Value this->maxInteger() means all documents.

    bool ok, has_grammar, uses_rle;

    const static usint RLE_FLAG = 0x01;

//--------------------------------------------------------------------------

    template<class T>
    std::vector<T>*
    documentListing(pair_type range) const
    {
      std::vector<T>* result = new std::vector<T>;
      DeltaVector::Iterator block_iter(*(this->leaf_ranges));

      // Process the part of the range before the first full block.
      pair_type first_block = block_iter.valueAfter(range.first);
      if(first_block.first > range.first)
      {
        usint temp = std::min(range.second, first_block.first - 1);
        this->processRange(pair_type(range.first, temp), *result);
        if(range.second <= temp) { return result; }
      }
      range.first = first_block.second;

      // Process the part of the range after the last full block.
      pair_type last_block = block_iter.valueBefore(range.second);
      pair_type next_block = block_iter.nextValue();
      if(range.second < next_block.first - 1)
      {
        this->processRange(pair_type(last_block.first, range.second), *result);
        if(first_block.second >= last_block.second) { return result; }
        last_block.second--;
      }
      range.second = last_block.second;

      // Process the full blocks.
      SuccinctVector::Iterator child_iter(*(this->first_children));
      usint prev_block = this->getNumberOfNodes(), run = 0;
      for(usint i = range.first; i <= range.second; )
      {
        usint next_leaf = i + 1;
        while(child_iter.isSet(i))  // Current block is the leftmost child of its parent.
        {
          pair_type parent = this->parentOf(i, child_iter);
          if(parent.second <= range.second + 1) // Range contains the subtree of the parent.
          {
            i = parent.first; next_leaf = parent.second;
          }
          else { break; }
        }
        if(i == prev_block + run) { run++; }
        else
        {
          if(this->addBlocks(prev_block, run, *result)) { delete result; return this->allDocuments<T>(); }
          prev_block = i; run = 1;
        }
        i = next_leaf;
      }
      if(this->addBlocks(prev_block, run, *result)) { delete result; return this->allDocuments<T>(); }

      return result;
    }

    template<class T>
    void
    processRange(pair_type range, std::vector<T>& result) const
    {
      if(isEmpty(range)) { return; }

      usint* res = this->rlcsa.locate(range);
      this->rlcsa.getSequenceForPosition(res, length(range));
      for(usint i = 0; i < length(range); i++)
      {
        this->addItem(res[i], result);
      }
    }

    template<class T>
    bool
    addBlocks(usint first, usint number, std::vector<T>& result) const
    {
      if(first >= this->getNumberOfNodes() || number == 0) { return false; }

      SuccinctVector::Iterator iter(*(this->block_borders));
      usint from = iter.select(first);
      usint to = (number == 1 ? iter.selectNext() : iter.select(first + number));

      for(usint i = from; i < to; i++)
      {
        usint val = this->blocks->readItemConst(i);
        if(val >= this->maxInteger()) { return true; }  // Block contains all documents.
        if(val < this->getNumberOfDocuments()) // Value is a document id.
        {
          this->addItem(val, result);
        }
        else  // Value is a rule id.
        {
          val -= this->getNumberOfDocuments();
          SuccinctVector::Iterator rule_iter(*(this->rule_borders));
          usint rfrom = rule_iter.select(val);
          usint rto = rule_iter.selectNext();

          if(this->usesRLE())
          {
            usint run_start = 0;
            bool want_run_length = false;
            for(usint j = rfrom; j < rto; j++)
            {
              val = this->rules->readItemConst(j);
              if(val >= this->getNumberOfDocuments()) { return true; }  // Rule covers all documents.
              else if(want_run_length) { this->addRun(run_start, val, result); want_run_length = false; }
              else { run_start = val; want_run_length = true; }
            }
          }
          else
          {
            for(usint j = rfrom; j < rto; j++)
            {
              val = this->rules->readItemConst(j);
              if(val >= this->getNumberOfDocuments()) { return true; }  // Rule covers all documents.
              else { this->addItem(val, result); }
            }
          }
        }
      }

      return false;
    }

    template<class T>
    std::vector<T>* allDocuments() const
    {
      std::vector<T>* result = new std::vector<T>;
      this->addRun(0, this->getNumberOfDocuments(), *result);
      return result;
    }

//--------------------------------------------------------------------------

    inline void addRun(usint from, usint length, std::vector<usint>& result) const
    {
      for(usint i = from; i < from + length; i++) { result.push_back(i); }
    }

    inline void addRun(usint from, usint length, std::vector<pair_type>& result) const
    {
      result.push_back(pair_type(from, from + length - 1));
    }

    inline void addItem(usint item, std::vector<usint>& result) const
    {
      result.push_back(item);
    }
    
    inline void addItem(usint item, std::vector<pair_type>& result) const
    {
      result.push_back(pair_type(item, item));
    }

//--------------------------------------------------------------------------

    inline pair_type parentOf(usint tree_node, SuccinctVector::Iterator& iter) const
    {
      usint par = this->parents->readItemConst(iter.rank(tree_node) - 1);
      usint next_leaf = this->next_leaves->readItemConst(par);
      return pair_type(par + this->getNumberOfLeaves(), next_leaf);
    }

//--------------------------------------------------------------------------

    // These are not allowed.
    DocArray();
    DocArray(const DocArray&);
    DocArray& operator = (const DocArray&);
};


} // namespace CSA


#endif // _RLCSA_DOCARRAY_H
