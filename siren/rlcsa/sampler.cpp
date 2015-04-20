#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif

#include "sampler.h"
#include "rlcsa.h"


namespace CSA
{

Sampler::Sampler(uint _size) :
  size(_size), status(NOT_READY),
  samples(0), number_of_samples(0)
{
}

Sampler::~Sampler()
{
  delete[] this->samples; this->samples = 0;
}

void
Sampler::writeTo(const std::string& base_name) const
{
  std::string sample_name = base_name + SA_SAMPLES_EXTENSION;
  std::ofstream sample_file(sample_name.c_str(), std::ios_base::binary);
  if(!sample_file)
  {
    std::cerr << "Sampler: Error creating sample file!" << std::endl;
    return;
  }

  for(uint i = 0; i < this->number_of_samples; i++)
  {
    sample_file.write((char*)&(this->samples[i].second), sizeof(uint));
  }
  sample_file.close();
}

pair_type*
Sampler::getSamples(short_pair* sa, uint number_of_sequences, usint threads)
{
  if(sa == 0 || this->samples == 0 || this->status != SAMPLED) { return 0; }

  for(uint i = 0; i < this->number_of_samples; i++)
  {
    this->samples[i].first = sa[this->samples[i].second].second - number_of_sequences;
  }
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif
  parallelSort(samples, samples + this->number_of_samples);

  pair_type* result = this->samples; this->samples = 0;
  this->status = NOT_READY;
  return result;
}

//--------------------------------------------------------------------------

WeightedSampler::WeightedSampler(weight_type* weights, uint _size, bool _use_psi) :
  Sampler(_size),
  use_psi(_use_psi), adjustment(0),
  path_weights(0), predecessors(0),
  edge_weights(0), edge_totals(0),
  nodes(0), node_positions(0)
{
  if(weights == 0 || this->size == 0)
  {
    delete[] weights; return;
  }

  this->calculateEdges(weights);
  this->path_weights = new sum_type[this->nodes];
  this->predecessors = new uint[this->nodes];

  this->status = READY;
}

WeightedSampler::~WeightedSampler()
{
  this->cleanUp();
  delete[] this->node_positions; this->node_positions = 0;
}

void
WeightedSampler::cleanUp()
{
  delete[] this->path_weights; this->path_weights = 0;
  delete[] this->predecessors; this->predecessors = 0;
  delete[] this->edge_weights; this->edge_weights = 0;
  delete[] this->edge_totals; this->edge_totals = 0;
}

//--------------------------------------------------------------------------

bool
WeightedSampler::buildSamples(uint sample_rate, sum_type initial_adjustment, usint threads)
{
  if(sample_rate == 0 || this->status != READY) { return false; }

  usint force = 0;
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  force = this->nodes / (threads * threads);
  #endif

  delete[] this->samples; this->samples = 0;
  this->number_of_samples = (this->size + sample_rate - 1) / sample_rate;

  sum_type average_weight;
  if(this->use_psi)
  {
    average_weight = this->edge_totals[this->nodes - 1] / this->size;
  }
  else
  {
    average_weight = this->edge_totals[0] / this->size;
  }
  average_weight = std::max(average_weight, (sum_type)1);

  // Should be -3 * max_weight, 3 * max_weight, if we allow negative weights.
  sum_type left = 1, right = std::numeric_limits<sum_type>::max();
  if(initial_adjustment > 0) { this->adjustment = initial_adjustment; }
  else { this->adjustment = average_weight; }
  uint shortest = 0, longest = 0;
  uint* shortest_path = 0;
  uint* longest_path = 0;

  // Find adjustment such that the shortest minimum weight path has at most
  // 'number' links, while the longest minimum weight path has at least 'number' links.
  // With adjustment 0, we can always find a path of weight 0: just take as many edges as
  // possible. Hence adjustment 0 will be used only when we can afford to sample every
  // suffix with a positive weight. This is equivalent to greedy sampling.
  std::cout << "left: " << left << ", right: " << right << std::endl;
  while(true)
  {
    std::cout << "Adjustment: " << this->adjustment << std::endl;

    shortest = this->minimumWeightPath(true, force);
    std::cout << "  shortest: " << shortest << " (weight " << this->path_weights[this->nodes - 1] << ")" << std::endl;
    if(shortest > this->number_of_samples)
    {
      left = this->adjustment + 1;
      this->adjustment = std::min(2 * this->adjustment, left + (right - left) / 2);
      continue;
    }
    shortest_path = new uint[shortest + 1];
    for(uint i = shortest + 1, cur = this->nodes - 1; i > 0; i--, cur = this->predecessors[cur])
    {
      shortest_path[i - 1] = cur;
    }

    longest = this->minimumWeightPath(false, force);
    std::cout << "  longest:  " << longest << " (weight " << this->path_weights[this->nodes - 1] << ")" << std::endl;
    if(longest < this->number_of_samples && this->adjustment > 1) // If adjustment will still be >= 1.
    {
      right = this->adjustment - 1;
      this->adjustment = left + (right - left) / 2;
      delete[] shortest_path; shortest_path = 0;
      continue;
    }
    longest_path = new uint[longest + 1];
    for(uint i = longest + 1, cur = this->nodes - 1; i > 0; i--, cur = this->predecessors[cur])
    {
      longest_path[i - 1] = cur;
    }

    break;
  }
  std::cout << std::endl;
  this->cleanUp();

  // Build a 'number'-link path from the saved paths.
  this->number_of_samples = std::min(this->number_of_samples, longest);
  this->buildPath(this->number_of_samples, longest_path, longest, shortest_path, shortest);
  delete[] shortest_path; shortest_path = 0;

  this->samples = new pair_type[this->number_of_samples];
  for(uint i = 0; i < this->number_of_samples; i++)
  {
    samples[i].second = this->node_positions[longest_path[i]];
  }
  delete[] longest_path; longest_path = 0;

  this->status = SAMPLED;
  return true;
}

//--------------------------------------------------------------------------

void
WeightedSampler::buildPath(uint links, uint* path_a, uint length_a, uint* path_b, uint length_b)
{
  uint cur_a = 0, cur_b = 0;
  uint adj = length_a - links;

  while(true)
  {
    if(cur_a == cur_b + adj && path_b[cur_b] <= path_a[cur_a] && path_b[cur_b + 1] >= path_a[cur_a + 1])
    {
      break;
    }

    if(path_a[cur_a] <= path_b[cur_b]) { cur_a++; }
    else                               { cur_b++; }
  }

  // Combine the prefix of 'path_b' until 'cur_b - 1' with the suffix of 'path_a' from 'cur_a'.
  for(uint i = 0; i <= cur_b; i++) { path_a[i] = path_b[i]; }
  for(uint i = cur_b + 1, j = cur_a + 1; j < length_a; i++, j++) { path_a[i] = path_a[j]; }
}

//--------------------------------------------------------------------------

uint
WeightedSampler::minimumWeightPath(bool shortest, uint force)
{
  if(force == 0) { force = this->nodes; }
  std::vector<uint> end_points;
  std::vector<pair_type> partial_results;
  for(uint i = 0; i < this->nodes - 1; i += force)
  {
    end_points.push_back(i);
  }
  end_points.push_back(this->nodes - 1);
  uint parts = end_points.size() - 1;
  for(uint i = 0; i < end_points.size(); i++)
  {
    this->path_weights[end_points[i]] = 0;
    this->predecessors[end_points[i]] = 0;
    partial_results.push_back(pair_type(0, 0));
  }

  #pragma omp parallel for schedule(static)
  for(uint part = 0; part < parts; part++)
  {
    std::deque<node_type> active_nodes;
    active_nodes.push_back(node_type(end_points[part], 0));

    for(uint i = end_points[part] + 1; i < end_points[part + 1]; i++)
    {
      node_type current = this->link(active_nodes, i, shortest);

      // Retire the active nodes that are before the linked predecessor.
      while(active_nodes.front().first < this->predecessors[i])
      {
        active_nodes.pop_front();
      }

      // Retire the first active nodes, if the following ones are strictly better.
      while(active_nodes.size() > 1 &&
            this->isStrictlyBetter(active_nodes[1].first, active_nodes[0].first, i + 1))
      {
        active_nodes.pop_front();
      }

      // Retire the last active nodes, if the preceding node or the current node is always
      // strictly better.
      while(active_nodes.size() > 1 &&
            this->bridge(active_nodes[active_nodes.size() - 2].first, active_nodes.back().first, i))
      {
        active_nodes.pop_back();
      }

      // Add the current node to active nodes, if it is at least as good as the active nodes
      // for some future node.
      if(this->isAsGood(i, active_nodes.back().first))
      {
        active_nodes.push_back(current);
      }
    }

    partial_results[part + 1] = this->finalLink(active_nodes, end_points[part + 1], shortest);
  }

  for(uint i = 1; i < partial_results.size(); i++)
  {
    partial_results[i].first += partial_results[i - 1].first;
    partial_results[i].second += partial_results[i - 1].second;
  }
  this->path_weights[this->nodes - 1] = partial_results[parts].second;

  return partial_results[parts].first;
}

node_type
WeightedSampler::link(std::deque<node_type>& active_nodes, uint to, bool shortest)
{
  this->path_weights[to] = std::numeric_limits<sum_type>::max();
  node_type current(to, 0);
  if(shortest) { current.second = std::numeric_limits<uint>::max(); }

  for(std::deque<node_type>::iterator iter = active_nodes.begin(); iter != active_nodes.end(); ++iter)
  {
    sum_type weight = this->getPathWeight((*iter).first, to);
    if(weight > this->path_weights[to]) { return current; }
    if(weight < this->path_weights[to] ||
       (shortest && (*iter).second + 1 < current.second) ||
       (!shortest && (*iter).second + 1 > current.second))
    {
      this->path_weights[to] = weight;
      this->predecessors[to] = (*iter).first;
      current.second = (*iter).second + 1;
    }
  }

  return current;
}

pair_type
WeightedSampler::finalLink(std::deque<node_type>& active_nodes, uint to, bool shortest)
{
  pair_type result(0, std::numeric_limits<sum_type>::max());
  result.first = (shortest ? std::numeric_limits<uint>::max() : 0);

  for(std::deque<node_type>::iterator iter = active_nodes.begin(); iter != active_nodes.end(); ++iter)
  {
    sum_type weight = this->getPathWeight((*iter).first, to);
    if(weight > result.second) { return result; }
    if(weight < result.second ||
       (shortest && (*iter).second + 1 < result.first) ||
       (!shortest && (*iter).second + 1 > result.first))
    {
      this->predecessors[to] = (*iter).first;
      result = pair_type((*iter).second + 1, weight);
    }
  }

  return result;
}

bool
WeightedSampler::bridge(uint a, uint b, uint c)
{
/*
  Note that 'c' is always less than 'this->size', and 'b' will not be added to active nodes,
  unless it offers a strictly better path to 'this->size' than 'a'.
*/

  uint low = c + 1, high = this->nodes - 1, mid = low;
  while(low < high)
  {
    if(this->isStrictlyBetter(a, b, mid))
    {
      low = mid + 1;
      mid = std::min(low + (high - low) / 2, c + 2 * (mid - c));
    }
    else { high = mid; mid = low + (high - low) / 2; }
  }
  // 'high' will now be the first position for which 'a' is not strictly better than 'b'.

  return this->isStrictlyBetter(c, b, high);
}

//--------------------------------------------------------------------------

void
WeightedSampler::calculateEdges(weight_type* weights)
{
  // Build nodes only for text positions with a positive weight.
  this->nodes = 2;
  for(uint i = 1; i < this->size; i++) { if(weights[i] > 0) { this->nodes++; } }
  this->node_positions = new uint[this->nodes];
  this->node_positions[0] = 0; this->node_positions[this->nodes - 1] = this->size;

  this->edge_totals  = new sum_type[this->nodes];
  this->edge_totals[0] = (this->use_psi ? 0 : weights[0]);
  for(uint i = 1, j = 1; i < this->size; i++)
  {
    if(weights[i] > 0) { this->edge_totals[j] = weights[i]; this->node_positions[j] = i; j++; }
  }
  this->edge_totals[this->nodes - 1] = 0;
  delete[] weights;
  this->edge_weights = new sum_type[this->nodes];


  if(this->use_psi)
  {
    this->edge_weights[0] = 0;
    for(uint i = 1; i < this->nodes; i++)
    {
      this->edge_weights[i] = this->edge_weights[i - 1] + this->getDistance(i - 1, i) * this->edge_totals[i - 1];
      this->edge_totals[i] += this->edge_totals[i - 1];
    }
  }
  else
  {
    this->edge_weights[this->nodes - 1] = 0;
    for(uint i = this->nodes - 1; i > 0; i--)
    {
      this->edge_weights[i - 1] = this->edge_weights[i] + this->getDistance(i - 1, i) * this->edge_totals[i];
      this->edge_totals[i - 1] += this->edge_totals[i];
    }
  }
}

sum_type
WeightedSampler::getEdgeWeight(uint from, uint to)
{
  if(this->use_psi)
  {
    return this->edge_weights[to] - this->edge_weights[from] -
           this->getDistance(from, to) * this->edge_totals[from] + this->adjustment;
  }
  else
  {
    return this->edge_weights[from] - this->edge_weights[to] -
           this->getDistance(from, to) * this->edge_totals[to] + this->adjustment;
  }
}

uint
WeightedSampler::getDistance(uint from, uint to)
{
  return this->node_positions[to] - this->node_positions[from];
}

sum_type
WeightedSampler::getPathWeight(uint through, uint to)
{
  return this->path_weights[through] + this->getEdgeWeight(through, to);
}

bool
WeightedSampler::isBetter(uint a, uint b, uint to)
{
  return (this->getPathWeight(a, to) <= this->getPathWeight(b, to));
}

bool
WeightedSampler::isStrictlyBetter(uint a, uint b, uint to)
{
  return (this->getPathWeight(a, to) < this->getPathWeight(b, to));
}

bool
WeightedSampler::isAsGood(uint a, uint b)
{
  sum_type a_val = this->path_weights[a], b_val = this->path_weights[b];

  if(this->use_psi)
  {
    a_val += this->getDistance(b, a) * this->edge_totals[b];
    b_val += this->edge_weights[a] - this->edge_weights[b];
    b_val += this->getDistance(a, this->nodes - 1) * (this->edge_totals[a] - this->edge_totals[b]);
  }
  else
  {
    a_val += this->getDistance(b, a) * this->edge_totals[this->nodes - 1];
    b_val += this->edge_weights[b] - this->edge_weights[a];
  }

  return (a_val <= b_val);
}

//--------------------------------------------------------------------------

SemiGreedySampler::SemiGreedySampler(weight_type* _weights, uint _size) :
  Sampler(_size),
  weights(0),
  sample_rate(0)
{
  if(_weights == 0 || this->size == 0)
  {
    delete[] weights; return;
  }

  this->weights = new std::pair<weight_type, uint>[this->size];
  for(uint i = 0; i < this->size; i++)
  {
    this->weights[i].first = _weights[i];
    this->weights[i].second = i;
  }
  delete[] _weights;
  this->status = READY;
}

SemiGreedySampler::~SemiGreedySampler()
{
  this->cleanUp();
}

void
SemiGreedySampler::cleanUp()
{
  delete[] this->weights; this->weights = 0;
}

bool
SemiGreedySampler::buildSamples(uint _sample_rate, double greediness)
{
  if(_sample_rate == 0 || this->status != READY || greediness < 0.0 || greediness > 1.0)
  {
    return false;
  }

  uint number = (this->size + _sample_rate - 1) / _sample_rate;
  this->samples = new pair_type[number];
  if(greediness == 1.0) { _sample_rate = this->size; }
  else
  {
    _sample_rate = std::min((double)(this->size), _sample_rate / (1.0 - greediness));
  }

  this->buildRegularSamples(_sample_rate);
  this->addGreedySamples(number);

  this->cleanUp();
  this->status = SAMPLED;
  return true;
}

void
SemiGreedySampler::buildRegularSamples(uint _sample_rate)
{
  this->sample_rate = _sample_rate;
  for(uint i = 0; i < this->size; i += this->sample_rate)
  {
    this->samples[this->number_of_samples].second = i;
    this->number_of_samples++;
  }
}

void
SemiGreedySampler::addGreedySamples(uint total_samples)
{
  if(this->number_of_samples >= total_samples) { return; }

  std::sort(this->weights, this->weights + this->size);
  for(uint i = this->size - 1; this->number_of_samples < total_samples; i--)
  {
    if(this->weights[i].second % this->sample_rate != 0)
    {
      if(this->weights[i].first == 0) { break; }  // The further samples will be weight 0.
      this->samples[this->number_of_samples].second = this->weights[i].second;
      this->number_of_samples++;
    }
  }
}

//--------------------------------------------------------------------------

}; // namespace CSA
