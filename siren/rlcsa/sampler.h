#ifndef SAMPLER_H
#define SAMPLER_H


#include <deque>

#include "misc/utils.h"


namespace CSA
{


typedef uint  weight_type;
typedef usint sum_type;

typedef std::pair<uint, uint> node_type;


class Sampler
{
  public:
    explicit Sampler(uint size);
    virtual ~Sampler();

    // Writes the sampled positions (as 32-bit uints) to the sample file.
    void writeTo(const std::string& base_name) const;

    // Returns the given number of samples as sorted (index, value) pairs.
    // The user must delete the samples.
    pair_type* getSamples(short_pair* sa, uint number_of_sequences, usint threads);

    const static uint NOT_READY = 0;
    const static uint READY = 1;
    const static uint SAMPLED = 2;

    inline uint getStatus() const { return this->status; }
    inline uint getSize() const { return this->size; }
    inline uint getItems() const { return this->number_of_samples; }

  protected:
    uint size, status;

    pair_type* samples;
    uint       number_of_samples;

    // These are not allowed.
    Sampler();
    Sampler(const Sampler&);
    Sampler& operator = (const Sampler&);
};


class WeightedSampler : public Sampler
{
  public:
    // We assume that weights are non-negative.
    // 'weights' will be deleted in the constructor.
    WeightedSampler(weight_type* weights, uint _size, bool _use_psi);
    virtual ~WeightedSampler();

    // Call this to build the samples. Returns 'true' if ok.
    // Set 'initial_adjustment' to nonzero if you know the correct value.
    bool buildSamples(uint sample_rate, sum_type initial_adjustment, usint threads);

  private:
    const static usint VECTOR_BLOCK_SIZE = 64;

    bool      use_psi;
    sum_type  adjustment;

    sum_type* path_weights;
    uint*     predecessors;

    /*
      When using LF: Edges (V_i, V_n)
        edge_weights[i] is \sum_{j = i}^{n - 1} distance(i, j) weights[j]
        edge_totals[i]  is \sum_{j = i}^{n - 1} weights[j]

      When using Psi: Edges (V_0, V_i)
        edge_weights[i] is \sum_{j = 1}^{i} distance(j, i) weights[j]
        edge_totals[i]  is \sum_{j = 1}^{i} weights[j]

      Note that edge_weights can overflow. We assume that sum_type is unsigned, so that
      modular arithmetic allows us to compute the weights of short enough edges, even if
      an overflow occurs.
    */
    sum_type* edge_weights;
    sum_type* edge_totals;  // FIXME do not store explicitly?

    uint  nodes;
    uint* node_positions;

    // Delete unnecessary structures to save space after sampling has been done.
    void cleanUp();

    // Build a path of 'links' edges that is no heavier than 'path_a' and 'path_b'.
    // Store it in 'path_a', which is assumed to be the longer one.
    void buildPath(uint links, uint* path_a, uint length_a, uint* path_b, uint length_b);

    // Finds the shortest or the longest path of minimum weight.
    // If we force one out of 'force' nodes to be included in the path, the minimum weight
    // path can be found in parallel.
    uint minimumWeightPath(bool shortest, uint force = 0);
    node_type link(std::deque<node_type>& active_nodes, uint to, bool shortest);
    pair_type finalLink(std::deque<node_type>& active_nodes, uint to, bool shortest);

    // Is 'a' or 'c' always strictly better than 'b'?
    bool bridge(uint a, uint b, uint c);

    void calculateEdges(weight_type* weights);

    sum_type getEdgeWeight(uint from, uint to);
    sum_type getPathWeight(uint through, uint to);

    uint getDistance(uint from, uint to);

    // Is the path through 'a' to 'to' (strictly) better than through 'b' to 'to'?
    bool isBetter(uint a, uint b, uint to);
    bool isStrictlyBetter(uint a, uint b, uint to);

    // Is the path through 'a' to end at least as good through 'b' to end?
    // We assume that a > b.
    // This is the only function that does not use getEdgeWeight() to get edge weights.
    // This allows us to use it for long edges.
    bool isAsGood(uint a, uint b);

    // These are not allowed.
    WeightedSampler();
    WeightedSampler(const WeightedSampler&);
    WeightedSampler& operator = (const WeightedSampler&);
};


class SemiGreedySampler : public Sampler
{
  public:
    // We assume that weights are non-negative.
    // 'weights' will be deleted in the constructor.
    SemiGreedySampler(weight_type* weights, uint _size);
    virtual ~SemiGreedySampler();

    // Call this to build the samples. Returns 'true' if ok.
    // Greediness should be from 0.0 to 1.0.
    bool buildSamples(uint _sample_rate, double greediness);

  private:
    std::pair<weight_type, uint>* weights;
    uint sample_rate;

    void cleanUp();

    void buildRegularSamples(uint _sample_rate);
    void addGreedySamples(uint total_samples);

    // These are not allowed.
    SemiGreedySampler();
    SemiGreedySampler(const SemiGreedySampler&);
    SemiGreedySampler& operator = (const SemiGreedySampler&);
};


}; // namespace CSA


#endif // SAMPLER_H
