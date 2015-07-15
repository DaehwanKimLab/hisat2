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

template <typename T, typename index_t>
struct Params {
	index_t     (*hash)(T&);
	T**         begin;
	int         log_size;
	int         num;
};

//basically used to wrap together calls to bin_sort
template <typename T, typename CMP, typename index_t>
void bin_sort_worker(void* vp) {
	Params<T, index_t>* params = (Params<T, index_t>*)vp;
	index_t (*hash)(T&) = params->hash;
	T**     begin       = params->begin;
	int     log_size    = params->log_size;
	int     num         = params->num;
	for(int i = 0; i < num; i++) {
		if(begin[i + 1] - begin[i] > 1) bin_sort<T, CMP, index_t>(begin[i], begin[i + 1], hash, log_size);
	}
}

template <typename T, typename CMP, typename index_t>
void bin_sort_no_copy(T* begin, T* end, index_t (*hash)(T&), index_t maxv, int nthreads = 1) {
	const int SHIFT = 7;
	const int BLOCKS = (1 << (SHIFT + 1));

	int log_size = sizeof(maxv) * 8;
	while(!((1 << log_size) & maxv)) log_size--;
	int right_shift = log_size - SHIFT;

	// {(maxv >> right_shift) + 1 <= BLOCKS},
	int occupied = (maxv >> right_shift) + 1;

	// count number in each bin
	index_t count[BLOCKS] = {0};
	for(T* curr = begin; curr != end; curr++) {
		count[hash(*curr) >> right_shift]++;
	}
	// sum numbers to create an index
	T* index[BLOCKS + 1];
	T* place[BLOCKS];
	index[0] = place[0] = begin;
	for(int i = 1; i < occupied; i++) {
		index[i] = place[i] = index[i - 1] + count[i - 1];
	}
	index[occupied] = end;
	//put objects in proper place
	for(int bin = 0; bin < occupied; bin++) {
		while(place[bin] != index[bin + 1]) {
			T curr = *place[bin];
			int x = hash(curr) >> right_shift;
			while(x != bin) { // switched inner loop here, removed branch statement
				T temp = *place[x];
				*place[x]++ = curr;
				curr = temp;
				x = hash(curr) >> right_shift;
			}
			*place[bin]++ = curr;
		}
	}
	//sort partitions
	if(nthreads == 1) {
		for(int bin = 0; bin < occupied; bin++) {
			if(index[bin + 1] - index[bin] > 1) bin_sort<T, CMP, index_t>(index[bin], index[bin + 1], hash, right_shift);
		}
	} else {
		AutoArray<tthread::thread*> threads(nthreads);
        EList<Params<T, index_t> > params; params.resizeExact(nthreads);
		int st = 0;
		for(int i = 0; i < nthreads; i++) {
			params[i].hash = hash;
			params[i].begin = index + st;
			params[i].log_size = right_shift;
			params[i].num = occupied / nthreads;
			threads[i] = new tthread::thread(&bin_sort_worker<T, CMP, index_t>, (void*)&params[i]);
			st += params[i].num;
		}
		//do any remaining bins using main thread
		for(int bin = st; bin < occupied; bin++) {
			if(index[bin + 1] - index[bin] > 1) bin_sort<T, CMP, index_t>(index[bin], index[bin + 1], hash, right_shift);
		}
		for(int i = 0; i < nthreads; i++) {
			threads[i]->join();
		}
	}
}

template <typename T, typename CMP, typename index_t>
void bin_sort_copy(T* begin, T* end, T* o, index_t (*hash)(T&), index_t maxv, int nthreads = 1) {
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
	if(nthreads == 1) {
		bin_sort<T, CMP, index_t>(o, index[0], hash, right_shift);
		for(int bin = 1; bin < occupied; bin++) {
			if(index[bin] - index[bin - 1] > 1) bin_sort<T, CMP, index_t>(index[bin - 1], index[bin], hash, right_shift);
		}
	} else {
		AutoArray<tthread::thread*> threads(nthreads);
        EList<Params<T, index_t> > params; params.resizeExact(nthreads);
		int st = 0;
		for(int i = 0; i < nthreads; i++) {
			params[i].hash = hash;
			params[i].begin = index + st;
			params[i].log_size = right_shift;
			params[i].num = occupied / nthreads;
			threads[i] = new tthread::thread(&bin_sort_worker<T, CMP, index_t>, (void*)&params[i]);
			st += params[i].num;
		}
		//do any remaining bins using main thread,
		for(int bin = st; bin < occupied - 1; bin++) {
			if(index[bin + 1] - index[bin] > 1) bin_sort<T, CMP, index_t>(index[bin], index[bin + 1], hash, right_shift);
		}
		bin_sort<T, CMP, index_t>(o, index[0], hash, right_shift);
		for(int i = 0; i < nthreads; i++) {
			threads[i]->join();
		}
	}
}

#endif //RADIX_SORT_H_
