adaptive_samples.o: adaptive_samples.cpp adaptive_samples.h rlcsa.h \
 bits/deltavector.h bits/bitvector.h bits/bitbuffer.h \
 bits/../misc/definitions.h bits/rlevector.h bits/nibblevector.h \
 bits/succinctvector.h sasamples.h sampler.h misc/utils.h \
 misc/definitions.h bits/bitbuffer.h alphabet.h misc/definitions.h \
 lcpsamples.h bits/array.h misc/parameters.h suffixarray.h
alphabet.o: alphabet.cpp alphabet.h misc/definitions.h
build_plcp.o: build_plcp.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
build_rlcsa.o: build_rlcsa.cpp rlcsa_builder.h rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
build_sa.o: build_sa.cpp suffixarray.h misc/definitions.h misc/utils.h \
 misc/definitions.h
display_test.o: display_test.cpp rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
docarray.o: docarray.cpp docarray.h rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
document_graph.o: document_graph.cpp rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h docarray.h
extract_sequence.o: extract_sequence.cpp rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
lcp_test.o: lcp_test.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
lcpsamples.o: lcpsamples.cpp lcpsamples.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/array.h misc/utils.h misc/definitions.h
locate_test.o: locate_test.cpp rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
parallel_build.o: parallel_build.cpp rlcsa_builder.h rlcsa.h \
 bits/deltavector.h bits/bitvector.h bits/bitbuffer.h \
 bits/../misc/definitions.h bits/rlevector.h bits/nibblevector.h \
 bits/succinctvector.h sasamples.h sampler.h misc/utils.h \
 misc/definitions.h bits/bitbuffer.h alphabet.h misc/definitions.h \
 lcpsamples.h bits/array.h misc/parameters.h suffixarray.h
read_bwt.o: read_bwt.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
rlcsa.o: rlcsa.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h bits/vectors.h
rlcsa_builder.o: rlcsa_builder.cpp rlcsa_builder.h rlcsa.h \
 bits/deltavector.h bits/bitvector.h bits/bitbuffer.h \
 bits/../misc/definitions.h bits/rlevector.h bits/nibblevector.h \
 bits/succinctvector.h sasamples.h sampler.h misc/utils.h \
 misc/definitions.h bits/bitbuffer.h alphabet.h misc/definitions.h \
 lcpsamples.h bits/array.h misc/parameters.h suffixarray.h
rlcsa_grep.o: rlcsa_grep.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
rlcsa_test.o: rlcsa_test.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h adaptive_samples.h docarray.h
sample_lcp.o: sample_lcp.cpp rlcsa.h bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/rlevector.h \
 bits/nibblevector.h bits/succinctvector.h sasamples.h sampler.h \
 misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
sampler.o: sampler.cpp sampler.h misc/utils.h misc/definitions.h rlcsa.h \
 bits/deltavector.h bits/bitvector.h bits/bitbuffer.h \
 bits/../misc/definitions.h bits/rlevector.h bits/nibblevector.h \
 bits/succinctvector.h sasamples.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
sampler_test.o: sampler_test.cpp rlcsa.h bits/deltavector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/rlevector.h bits/nibblevector.h bits/succinctvector.h sasamples.h \
 sampler.h misc/utils.h misc/definitions.h bits/bitbuffer.h alphabet.h \
 misc/definitions.h lcpsamples.h bits/array.h misc/parameters.h \
 suffixarray.h
sasamples.o: sasamples.cpp sasamples.h sampler.h misc/utils.h \
 misc/definitions.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/deltavector.h bits/bitvector.h bits/bitbuffer.h
ss_test.o: ss_test.cpp misc/utils.h misc/definitions.h
suffixarray.o: suffixarray.cpp misc/utils.h misc/definitions.h \
 suffixarray.h misc/definitions.h
array.o: bits/array.cpp bits/array.h bits/../misc/definitions.h \
 bits/bitbuffer.h
bitbuffer.o: bits/bitbuffer.cpp bits/bitbuffer.h \
 bits/../misc/definitions.h
bitvector.o: bits/bitvector.cpp bits/bitvector.h bits/bitbuffer.h \
 bits/../misc/definitions.h
deltavector.o: bits/deltavector.cpp bits/deltavector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h
multiarray.o: bits/multiarray.cpp bits/multiarray.h bits/array.h \
 bits/../misc/definitions.h bits/bitbuffer.h bits/succinctvector.h \
 bits/bitvector.h
nibblevector.o: bits/nibblevector.cpp bits/nibblevector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/../misc/utils.h bits/../misc/definitions.h
rlevector.o: bits/rlevector.cpp bits/rlevector.h bits/bitvector.h \
 bits/bitbuffer.h bits/../misc/definitions.h bits/../misc/utils.h \
 bits/../misc/definitions.h
succinctvector.o: bits/succinctvector.cpp bits/succinctvector.h \
 bits/bitvector.h bits/bitbuffer.h bits/../misc/definitions.h \
 bits/../misc/utils.h bits/../misc/definitions.h
parameters.o: misc/parameters.cpp misc/parameters.h misc/definitions.h
utils.o: misc/utils.cpp misc/utils.h misc/definitions.h
convert_patterns.o: utils/convert_patterns.cpp \
 utils/../misc/definitions.h utils/../misc/utils.h \
 utils/../misc/definitions.h
extract_text.o: utils/extract_text.cpp utils/../misc/definitions.h
sort_wikipedia.o: utils/sort_wikipedia.cpp utils/../misc/utils.h \
 utils/../misc/definitions.h
split_text.o: utils/split_text.cpp utils/../misc/definitions.h \
 utils/../misc/utils.h utils/../misc/definitions.h
