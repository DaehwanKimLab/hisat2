#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "docarray.h"
#include "misc/utils.h"


namespace CSA
{

//--------------------------------------------------------------------------

STNode::STNode(uint lcp, pair_type sa_range) :
  string_depth(lcp), range(sa_range),
  parent(0), child(0), sibling(0), next(0),
  stored_documents(0), id(0), contains_all(false),
  docs(0)
{
}

STNode::~STNode()
{
  this->parent = 0;
  delete this->docs; this->docs = 0;
  delete this->child; this->child = 0;
  while(this->sibling != 0)
  {
    STNode* temp = this->sibling; this->sibling = this->sibling->sibling;
    temp->sibling = 0; delete temp;
  }
}

void
STNode::addChild(STNode* node)
{
  node->parent = this;
  if(this->child == 0) { this->child = node; return; }
  this->child->addSibling(node);
}

void
STNode::addSibling(STNode* node)
{
  if(this->sibling == 0) { this->sibling = node; return; }
  this->sibling->addSibling(node);
}

void
STNode::deleteChildren()
{
  delete this->child; this->child = 0;
}

void
STNode::addLeaves()
{
  std::stack<STNode*> nodestack;
  nodestack.push(this);

  while(!(nodestack.empty()))
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    if(curr->child == 0) { continue; }
    usint expect = curr->range.first;
    STNode* prev = 0;
    for(STNode* temp = curr->child; temp != 0; temp = temp->sibling)
    {
      while(expect < temp->range.first) { prev = curr->addLeaf(prev, temp, expect); expect++; }
      nodestack.push(temp);
      expect = temp->range.second + 1;
      prev = temp;
    }
    while(expect < curr->range.second + 1) { prev = curr->addLeaf(prev, 0, expect); expect++; }
  }
}

STNode*
STNode::addLeaf(STNode* left, STNode* right, usint pos)
{
  STNode* leaf = new STNode(0, pair_type(pos, pos));
  leaf->parent = this;
  if(left == 0) { this->child = leaf; }
  else          { left->sibling = leaf; }
  leaf->sibling = right;
  return leaf;
}

bool
STNode::verifyTree()
{
  std::stack<STNode*> nodestack;
  nodestack.push(this);

  while(!(nodestack.empty()))
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    if(curr->child == 0) { continue; }
    usint expect = curr->range.first, last = 0;
    for(STNode* temp = curr->child; temp != 0; temp = temp->sibling)
    {
      if(temp->range.first != expect)
      {
        std::cerr << "Error: node " << *curr << ", expected " << expect << ", got " << *temp << std::endl;
        std::cerr << "Children:";
        for(temp = curr->child; temp != 0; temp = temp->sibling) { std::cerr << " " << *temp; }
        std::cerr << std::endl;
        if(curr->parent) { std::cerr << "Parent: " << *(curr->parent) << std::endl; }
        return false;
      }
      nodestack.push(temp);
      expect = temp->range.second + 1;
      last = temp->range.second;
    }
    if(last != curr->range.second)
    {
      std::cerr << "Error: node " << *curr << ", children end at " << last << std::endl;
      return false;
    }
  }

  return true;
}

void
STNode::remove()
{
  if(this->child == 0 || this->parent == 0) { return; }  // Cannot remove root or leaves.

  STNode* curr;
  for(curr = this->child; curr != 0; curr = curr->sibling) { curr->parent = this->parent; }
  if(this == this->parent->child) { this->parent->child = this->child; }
  else
  {
    for(curr = this->parent->child; curr->sibling != this; curr = curr->sibling);
    curr->sibling = this->child;
  }
  for(curr = this->child; curr->sibling != 0; curr = curr->sibling);
  curr->sibling = this->sibling;

  this->parent = this->child = this->sibling = 0;
}

void
STNode::determineSize(uint& nodes, uint& leaves)
{
  nodes = 0; leaves = 0;
  std::stack<STNode*> nodestack;
  nodestack.push(this);

  while(!(nodestack.empty()))
  {
    STNode* curr = nodestack.top(); nodestack.pop(); nodes++;
    if(curr->child == 0) { leaves++; }
    for(STNode* temp = curr->child; temp != 0; temp = temp->sibling) { nodestack.push(temp); }
  }
}

void
STNode::computeStoredDocuments(usint* documents)
{
  this->stored_documents = 0;
  delete this->docs; this->docs = 0;
  if(this->child == 0)
  {
    this->docs = new std::vector<uint>(documents + this->range.first, documents + this->range.second + 1);
  }
  else
  {
    this->docs = new std::vector<uint>;
    for(STNode* curr = this->child; curr != 0; curr = curr->sibling)
    {
      this->docs->insert(this->docs->end(), curr->docs->begin(), curr->docs->end());
      this->stored_documents += curr->stored_documents;
    }
  }
  removeDuplicates(this->docs, (this->docs->size() >= DocArray::PARALLEL_SORT_THRESHOLD));
}

void
STNode::containsAllDocuments()
{
  if(this->docs != 0)
  {
    this->stored_documents = this->docs->size();
    delete this->docs; this->docs = 0;
  }
  for(STNode* curr = this; curr != 0; curr = curr->parent)
  {
    curr->contains_all = true;
  }
}

void
STNode::storeThisSet()
{
  if(!(this->contains_all)) { this->stored_documents = this->docs->size(); }
  for(STNode* curr = this->child; curr != 0; curr = curr->sibling)
  {
    delete curr->docs; curr->docs = 0;
  }
}

void
STNode::setNext()
{
  std::stack<STNode*> nodestack;
  std::vector<STNode*> previous;
  for(STNode* temp = this; temp != 0; temp = temp->child) { nodestack.push(temp); }
  this->next = 0;

  while(nodestack.top() != this)
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    for(STNode* temp = curr->sibling; temp != 0; temp = temp->child) { nodestack.push(temp); }

    if(curr->child == 0)
    {
      for(std::vector<STNode*>::iterator iter = previous.begin(); iter != previous.end(); ++iter)
      {
        (*iter)->next = curr;
      }
      previous.clear();
    }
    curr->next = 0; previous.push_back(curr);
  }
}

std::ostream&
operator<<(std::ostream& stream, const STNode& data)
{
  stream << "(depth = " << data.string_depth << ", range = " << data.range << ", parent = " << data.parent;

  stream << ", siblings = (";
  for(STNode* curr = data.sibling; curr != 0; curr = curr->sibling)
  {
    if(curr != data.sibling) { stream << ", "; }
    stream << curr;
  }
  stream << ")";

  stream << ", children = (";
  for(STNode* curr = data.child; curr != 0; curr = curr->sibling)
  {
    if(curr != data.child) { stream << ", "; }
    stream << curr;
  }
  stream << ")";

  stream << ")";
  return stream;
}

//--------------------------------------------------------------------------

DocArray::DocArray(const RLCSA& _rlcsa) :
  rlcsa(_rlcsa),
  leaf_ranges(0), first_children(0), parents(0), next_leaves(0),
  rule_borders(0), rules(0),
  block_borders(0), blocks(0),
  ok(false), has_grammar(false), uses_rle(false)
{
}

DocArray::DocArray(STNode* root, const RLCSA& _rlcsa) :
  rlcsa(_rlcsa),
  leaf_ranges(0), first_children(0), parents(0), next_leaves(0),
  rule_borders(0), rules(0),
  block_borders(0), blocks(0),
  ok(false), has_grammar(false), uses_rle(false)
{
  if(root == 0) { return; }

  uint nodes = 0, leaves = 0;
  root->determineSize(nodes, leaves);
  root->setNext();

  DeltaVector::Encoder range_encoder(RANGE_BLOCK_SIZE);
  SuccinctVector::Encoder child_encoder(PARENT_BLOCK_SIZE,
    nextMultipleOf(PARENT_BLOCK_SIZE, BITS_TO_BYTES(nodes)));

  std::stack<STNode*> nodestack;
  for(STNode* curr = root; curr != 0; curr = curr->child) { nodestack.push(curr); }
  while(nodestack.top() != root)
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    for(STNode* temp = curr->sibling; temp != 0; temp = temp->child) { nodestack.push(temp); }
    if(curr->child == 0) { range_encoder.addBit(curr->range.first); }
    if(curr == curr->parent->child) { child_encoder.setBit(curr->id); }
  }
  this->leaf_ranges = new DeltaVector(range_encoder, root->range.second + 1);
  this->first_children = new SuccinctVector(child_encoder, nodes);

  // Pointers to the parents for first children and to the next leaf for internal nodes.
  WriteBuffer par_buffer(nodes - leaves, length(nodes - leaves - 1));
  WriteBuffer next_buffer(nodes - leaves, length(leaves));
  SuccinctVector::Iterator iter(*(this->first_children));
  while(!(nodestack.empty())) { nodestack.pop(); }
  for(STNode* curr = root; curr != 0; curr = curr->child) { nodestack.push(curr); }
  while(nodestack.top() != root)
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    for(STNode* temp = curr->sibling; temp != 0; temp = temp->child) { nodestack.push(temp); }
    if(curr == curr->parent->child)
    {
      par_buffer.goToItem(iter.rank(curr->id) - 1);
      par_buffer.writeItem(curr->parent->id - leaves);
    }
    if(curr->id >= leaves)
    {
      next_buffer.writeItem(curr->next == 0 ? this->getNumberOfLeaves() : curr->next->id);
    }
  }
  this->parents = par_buffer.getReadBuffer();
  this->next_leaves = next_buffer.getReadBuffer();

  this->ok = true;
}

DocArray::DocArray(const RLCSA& _rlcsa, const std::string& base_name, bool load_grammar) :
  rlcsa(_rlcsa),
  leaf_ranges(0), first_children(0), parents(0), next_leaves(0),
  rule_borders(0), rules(0),
  block_borders(0), blocks(0),
  ok(false), has_grammar(false), uses_rle(false)
{
  std::string input_name = base_name + DOCUMENT_EXTENSION;
  std::ifstream input(input_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "DocArray: Error opening input file " << input_name << std::endl;
    return;
  }

  usint flags = 0;
  input.read((char*)(&flags), sizeof(flags));
  if(flags & RLE_FLAG) { this->uses_rle = true; }

  this->leaf_ranges = new DeltaVector(input);
  this->first_children = new SuccinctVector(input);

  this->parents = new ReadBuffer(input, this->getNumberOfInternalNodes(),
    length(this->getNumberOfInternalNodes() - 1));
  this->next_leaves = new ReadBuffer(input, this->getNumberOfInternalNodes(),
    length(this->getNumberOfLeaves()));
  this->ok = true;

  if(load_grammar)
  {
    this->rule_borders = new SuccinctVector(input);
    this->rules = new ReadBuffer(input, this->rule_borders->getSize(),
      length(this->getNumberOfDocuments()));
    this->block_borders = new SuccinctVector(input);
    this->blocks = new ReadBuffer(input, this->block_borders->getSize(),
      length(this->maxInteger()));
    this->has_grammar = true;
  }

  input.close();
}

DocArray::~DocArray()
{
  delete this->leaf_ranges; this->leaf_ranges = 0;
  delete this->first_children; this->first_children = 0;
  delete this->parents; this->parents = 0;
  delete this->next_leaves; this->next_leaves = 0;

  delete this->rule_borders; this->rule_borders = 0;
  delete this->rules; this->rules = 0;

  delete this->block_borders; this->block_borders = 0;
  delete this->blocks; this->blocks = 0;
}

//--------------------------------------------------------------------------

void
parseLine(std::string& row, std::vector<usint>& before, int separator, std::vector<usint>& after)
{
  bool found = false, num = false;
  usint value = 0;
  for(usint i = 0; i < row.length(); i++)
  {
    if(isdigit(row[i])) { value = 10 * value + row[i] - '0'; num = true; }
    else if(num)
    {
      if(found) { after.push_back(value); }
      else      { before.push_back(value); }
      value = 0; num = false;
    }
    if(row[i] == separator) { found = true; }
  }
  if(num)
  {
    if(found) { after.push_back(value); }
    else      { before.push_back(value); }
    value = 0; num = false;
  }
}

usint
parseSingletonFile(const std::string& file_name, const DocArray& docarray, std::vector<pair_type>& occurrences, bool print)
{
  std::ifstream input(file_name.c_str(), std::ios_base::binary);
  if(!input) { return 0; }

  usint total_singletons = 0;
  std::vector<std::string> rows;
  readRows(input, rows, true);
  input.close();

  bool invalid_blocks = false, invalid_docs = false;
  for(std::vector<std::string>::iterator iter = rows.begin(); iter != rows.end(); ++iter)
  {
    std::vector<usint> before, after;
    parseLine(*iter, before, ':', after);

    if(before.size() != 1 || !(docarray.nodeIsBlock(before[0])))
    {
      invalid_blocks = true;
      continue;
    }
    usint block_id = docarray.nodeToBlock(before[0]);

    usint prev_singleton = 0;
    for(std::vector<usint>::iterator iter = after.begin(); iter != after.end(); ++iter)
    {
      if(docarray.nodeIsDoc(*iter))
      {
        // There may be duplicates due to the implementation used to find bicliques.
        if(*iter != prev_singleton)
        {
          occurrences.push_back(pair_type(block_id, docarray.nodeToDoc(*iter)));
          total_singletons++; prev_singleton = *iter;
        }
      }
      else { invalid_docs = true; }
    }
  }

  if(print)
  {
    std::cout << "Singleton file" << file_name << ": " << total_singletons << " singletons" << std::endl;
    if(invalid_blocks)     { std::cout << "There were invalid block ids." << std::endl; }
    if(invalid_docs)       { std::cout << "There were invalid document ids." << std::endl; }
    std::cout << std::endl;
  }

  return total_singletons;
}

void
DocArray::readRules(const std::string& name_prefix, bool print)
{
  std::vector<pair_type> occurrences; // (block_id, rule_id)
  std::vector<usint> docs;
  SuccinctVector::Encoder rb_encoder(RULE_BLOCK_SIZE);
  SuccinctVector::Encoder rb_encoder_rle(RULE_BLOCK_SIZE);
  usint files = 0, current_rule = this->getNumberOfDocuments();
  usint total_rules = 0, total_occs = 0, total_size = 0, total_runs = 0;


  // Parse the individual files containing rules / bicliques.
  while(true)
  {
    std::ostringstream fn;
    fn << name_prefix << "-biclique-it-" << files << ".txt";
    std::ifstream input(fn.str().c_str(), std::ios_base::binary);
    if(!input) { break; }

    std::vector<std::string> rows;
    readRows(input, rows, true);
    input.close();

    bool invalid_blocks = false, nonoccurring_rules = false, invalid_docs = false, empty_rules = false;
    usint curr_rules = 0, curr_occs = 0, curr_size = 0, curr_runs = 0;
    for(std::vector<std::string>::iterator siter = rows.begin(); siter != rows.end(); ++siter)
    {
      std::vector<usint> before, after;
      parseLine(*siter, before, '-', after);

      // Block ids.
      bool occurs = false;
      for(std::vector<usint>::iterator iter = before.begin(); iter != before.end(); ++iter)
      {
        if(this->nodeIsBlock(*iter))
        {
          occurrences.push_back(pair_type(this->nodeToBlock(*iter), current_rule));
          curr_occs++; occurs = true;
        }
        else { invalid_blocks = true; }
      }
      if(!occurs) { nonoccurring_rules = true; }
      curr_rules++; current_rule++;

      // Document ids.
      rb_encoder.addBit(total_size + curr_size);
      rb_encoder_rle.addBit(docs.size());
      occurs = false;
      usint prev_doc = WORD_MAX - 1, run_length = 0;
      for(std::vector<usint>::iterator iter = after.begin(); iter != after.end(); ++iter)
      {
        if(this->nodeIsDoc(*iter))
        {
          usint doc_id = this->nodeToDoc(*iter);
          if(doc_id == prev_doc + 1) { run_length++; prev_doc++; }
          else
          {
            if(run_length > 0) { docs.push_back(run_length); }
            docs.push_back(doc_id); run_length = 1; curr_runs++;
            prev_doc = doc_id;
          }
          curr_size++; occurs = true;
        }
        else { invalid_docs = true; }
      }
      if(!occurs)
      {
        docs.push_back(this->getNumberOfDocuments());
        curr_size++; curr_runs++; empty_rules = true;
      }
      else { docs.push_back(run_length); }
    }

    if(print)
    {
      std::cout << "File " << files << ": "
        << curr_rules << " rules, "
        << curr_occs << " occurrences, "
        << curr_size << " document ids in "
        << curr_runs << " runs" << std::endl;
      if(invalid_blocks)     { std::cout << "There were invalid block ids." << std::endl; }
      if(nonoccurring_rules) { std::cout << "There were non-occurring rules." << std::endl; }
      if(invalid_docs)       { std::cout << "There were invalid document ids." << std::endl; }
      if(empty_rules)        { std::cout << "There were empty rules." << std::endl; }
      std::cout << std::endl;
    }
    files++;
    total_rules += curr_rules; total_occs += curr_occs; total_size += curr_size; total_runs += curr_runs;
  }


  // Parse the singleton file and the remaining part of the graph.
  usint total_singletons = 0, temp = 0;

  std::ostringstream ffn;
  ffn << name_prefix << "-it-" << (files - 1);
  std::string remainder_name = ffn.str();
  temp = parseSingletonFile(remainder_name, *this, occurrences, print);
  if(temp > 0) { total_singletons += temp; files++; }

  std::string singleton_name = name_prefix + ".singletons";
  temp = parseSingletonFile(singleton_name, *this, occurrences, print);
  if(temp > 0) { total_singletons += temp; files++; }


  // Build the structures.
  if(total_size < docs.size())
  {
    this->uses_rle = false;
    this->rule_borders = new SuccinctVector(rb_encoder, total_size);
    WriteBuffer rbuf(total_size, length(this->getNumberOfDocuments()));
    bool need_run_length = false; usint prev = 0;
    for(std::vector<usint>::iterator iter = docs.begin(); iter != docs.end(); ++iter)
    {
      if(!need_run_length)
      {
        prev = *iter;
        rbuf.writeItem(*iter);
        need_run_length = (*iter != this->getNumberOfDocuments());
      }
      else
      {
        for(usint i = 1; i < *iter; i++) { rbuf.writeItem(prev + i); }
        need_run_length = false;
      }
    }
    this->rules = rbuf.getReadBuffer();
  }
  else
  {
    this->uses_rle = true;
    this->rule_borders = new SuccinctVector(rb_encoder_rle, docs.size());
    WriteBuffer rbuf(docs.size(), length(this->getNumberOfDocuments()));
    for(std::vector<usint>::iterator iter = docs.begin(); iter != docs.end(); ++iter)
    {
      rbuf.writeItem(*iter);
    }
    this->rules = rbuf.getReadBuffer();
  }

  parallelSort(occurrences.begin(), occurrences.end());
  usint real_occs = occurrences.size();
  for(usint i = 0, looking_for = 0; i < occurrences.size(); i++)
  {
    while(looking_for < occurrences[i].first) { real_occs++; looking_for++; }
    looking_for = occurrences[i].first + 1;
  }
  real_occs += this->getNumberOfNodes() - 1 - occurrences.rbegin()->first;

  SuccinctVector::Encoder bb_encoder(BLOCK_BLOCK_SIZE);
  WriteBuffer bbuf(real_occs, length(this->maxInteger()));
  usint pos = 0;
  for(usint i = 0, looking_for = 0; i < occurrences.size(); i++)
  {
    while(looking_for < occurrences[i].first)
    {
      bb_encoder.addBit(pos);
      bbuf.writeItem(this->maxInteger()); pos++;
      looking_for++;
    }
    if(i == 0 || occurrences[i].first != occurrences[i - 1].first)
    {
      bb_encoder.addBit(pos);
      looking_for++;
    }
    bbuf.writeItem(occurrences[i].second); pos++;
  }
  for(usint i = occurrences.rbegin()->first + 1; i < this->getNumberOfNodes(); i++)
  {
    bb_encoder.addBit(pos);
    bbuf.writeItem(this->maxInteger()); pos++;
  }
  this->block_borders = new SuccinctVector(bb_encoder, real_occs);
  this->blocks = bbuf.getReadBuffer();

  if(print)
  {
    std::cout << files << " files: "
      << total_rules << " rules, "
      << total_occs << " occurrences, "
      << total_size << " document ids in "
      << total_runs << " runs, "
      << total_singletons << " singletons, "
      << (real_occs - total_occs - total_singletons) << " complete blocks" << std::endl;
    if(this->usesRLE()) { std::cout << "Using run-length encoded rules." << std::endl; }
    std::cout << std::endl;
  }
  this->has_grammar = true;
}

//--------------------------------------------------------------------------

void
DocArray::writeTo(const std::string& base_name) const
{
  std::string output_name = base_name + DOCUMENT_EXTENSION;
  std::ofstream output(output_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "DocArray: Error creating output file " << output_name << std::endl;
    return;
  }

  usint flags = 0; if(this->usesRLE()) { flags |= RLE_FLAG; }
  output.write((char*)(&flags), sizeof(flags));

  this->leaf_ranges->writeTo(output);
  this->first_children->writeTo(output);
  this->parents->writeBuffer(output);
  this->next_leaves->writeBuffer(output);

  if(this->hasGrammar())
  {
    this->rule_borders->writeTo(output);
    this->rules->writeBuffer(output);
    this->block_borders->writeTo(output);
    this->blocks->writeBuffer(output);
  }

  output.close();
}

usint
DocArray::reportSize(bool print) const
{
  usint bytes = sizeof(*this);

  usint tree_bytes = this->leaf_ranges->reportSize() +
    this->first_children->reportSize() +
    this->parents->reportSize() +
    this->next_leaves->reportSize();
  if(print)
  {
    std::cout << "Tree:            " << (tree_bytes / (double)MEGABYTE) << " MB" << std::endl;
  }

  usint rule_bytes = 0, block_bytes = 0;
  if(this->hasGrammar())
  {
    rule_bytes = this->rule_borders->reportSize() + this->rules->reportSize();
    block_bytes = this->block_borders->reportSize() + this->blocks->reportSize();
  }
  if(print && this->hasGrammar())
  {
    std::cout << "Rules:           " << (rule_bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Blocks:          " << (block_bytes / (double)MEGABYTE) << " MB" << std::endl;
  }

  bytes += tree_bytes + rule_bytes + block_bytes;
  if(print)
  {
    std::cout << "Total size:      " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }

  return bytes;
}

//--------------------------------------------------------------------------

std::vector<usint>*
DocArray::listDocuments(const std::string& pattern) const
{
  return this->listDocuments(this->rlcsa.count(pattern));
}

std::vector<usint>*
DocArray::listDocuments(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->getSize()) { return 0; }

  std::vector<usint>* result = this->documentListing<usint>(range);
  if(result != 0) { removeDuplicates(result, false); }

  return result;
}

usint
DocArray::count(const std::string& pattern) const
{
  if(!(this->isOk())) { return 0; }
  return this->count(this->rlcsa.count(pattern));
}

usint
DocArray::count(pair_type sa_range) const
{
  if(!(this->isOk()) || CSA::isEmpty(sa_range) || sa_range.second >= this->rlcsa.getSize()) { return 0; }

  std::vector<usint>* result = this->documentListing<usint>(sa_range);
  if(result == 0) { return 0; }
  removeDuplicates(result, false);
  usint num = result->size();
  delete result; result = 0;

  return num;
}

std::vector<pair_type>*
DocArray::listDocumentsRLE(const std::string& pattern) const
{
  return this->listDocumentsRLE(this->rlcsa.count(pattern));
}

std::vector<pair_type>*
DocArray::listDocumentsRLE(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->getSize()) { return 0; }

  std::vector<pair_type>* result = this->documentListing<pair_type>(range);
  if(result != 0) { mergeRanges(result, false); }

  return result;
}

std::vector<usint>*
DocArray::directListing(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->getSize()) { return 0; }
  std::vector<usint>* result = new std::vector<usint>;
  this->processRange(range, *result);
  removeDuplicates(result, false);
  return result;
}

std::vector<pair_type>*
DocArray::directListingRLE(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->getSize()) { return 0; }
  std::vector<pair_type>* result = new std::vector<pair_type>;
  this->processRange(range, *result);
  mergeRanges(result, false);
  return result;
}

//--------------------------------------------------------------------------

} // namespace CSA
