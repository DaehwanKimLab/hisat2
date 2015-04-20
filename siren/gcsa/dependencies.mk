build_automaton.o: build_automaton.cpp graph.h
build_index.o: build_index.cpp gcsa.h graph.h
clean_alignment.o: clean_alignment.cpp
determinize.o: determinize.cpp graph.h
gcsa.o: gcsa.cpp gcsa.h graph.h
gcsa_test.o: gcsa_test.cpp gcsa.h graph.h bwasearch.h parameter_handler.h \
  pattern_classifier.h
generator.o: generator.cpp
graph.o: graph.cpp graph.h gcsa.h
parameter_handler.o: parameter_handler.cpp parameter_handler.h
rlcsa_test.o: rlcsa_test.cpp bwasearch.h parameter_handler.h \
  pattern_classifier.h
transpose.o: transpose.cpp
verify.o: verify.cpp gcsa.h graph.h
