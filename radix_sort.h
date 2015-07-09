#ifndef RADIX_SORT_H_
#define RADIX_SORT_H_

template <typename T, typename CMP, typename index_t>
void bin_sort(T* begin, T* end, index_t (*hash)(T&), int log_size) {
	const int SHIFT = 7;
	const int BLOCKS = (1 << (SHIFT + 1));
	const int BLOCKS_MASK = BLOCKS - 1;

	if(end - begin < 2000 || !log_size) { //picked arbitrarily
		if(end > begin + 1) sort(begin, end, CMP());
		return;
	}

	// compute maximum of log_size - 7 and 0
	int right_shift = (log_size - SHIFT) * (log_size > SHIFT);
	// count number in each bin
	index_t count[BLOCKS] = {0};
	for(T* curr = begin; curr != end; curr++) {
		count[(hash(*curr) >> right_shift) & BLOCKS_MASK]++;
	}
	// sum numbers to create an index
	T* index[BLOCKS + 1];
	T* place[BLOCKS];
	index[0] = place[0] = begin;
	for(int i = 1; i < BLOCKS; i++) {
		index[i] = place[i] = index[i - 1] + count[i - 1];
	}
	index[BLOCKS] = end;
	//put objects in proper place
	for(int bin = 0; bin < BLOCKS; bin++) {
		while(place[bin] != index[bin + 1]) {
			T curr = *place[bin];
			int x = (hash(curr) >> right_shift) & BLOCKS_MASK;
			while(x != bin) { // switched inner loop here, removed branch statement
				T temp = *place[x];
				*place[x]++ = curr;
				curr = temp;
				x = (hash(curr) >> right_shift) & BLOCKS_MASK;
			}
			*place[bin]++ = curr;
		}
	}
	//sort partitions
	for(int bin = 0; bin < BLOCKS; bin++) {
		if(index[bin + 1] - index[bin] > 1) bin_sort<T, CMP, index_t>(index[bin], index[bin + 1], hash, right_shift);
	}
}

template <typename T, typename CMP, typename index_t>
void bin_sort_copy(T* begin, T* end, T* o, index_t (*hash)(T&), index_t maxv) {
	const int SHIFT = 7;
	const int BLOCKS = (1 << (SHIFT + 1));

	int log_size = sizeof(maxv) * 8;
	while(!((1 << log_size) & maxv)) log_size--;
	int right_shift = log_size - SHIFT;
	int occupied = (maxv >> right_shift) + 1;
	index_t count[BLOCKS] = {0};
	// count number in each bin ~ 1/5 of time
	for(T* curr = begin; curr != end; curr++) {
		count[hash(*curr) >> right_shift]++;
	}
	// sum numbers to create an index ~ takes very little time
	T* index[BLOCKS];
	index[0] = o;
	for(int i = 1; i < occupied; i++) {
		index[i] = index[i - 1] + count[i - 1];
	}
	// hash objects ~ 4/5 of time
	for(T* curr = begin; curr != end; curr++) {
		*index[hash(*curr) >> right_shift]++ = *curr;
	}
	//sort partitions
	bin_sort<T, CMP, index_t>(o, index[0], hash, right_shift);
	for(int bin = 1; bin < occupied; bin++) {
		if(index[bin] - index[bin - 1] > 1) bin_sort<T, CMP, index_t>(index[bin - 1], index[bin], hash, right_shift);
	}
}

#endif //RADIX_SORT_H_
