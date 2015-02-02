CC = g++
CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11
CPPFLAGS = -I$(TBBROOT)/include
# this is very unportable
LDFLAGS := $(LDFLAGS) $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib
LDLIBS = -ltbb -lboost_timer -lboost_chrono -lboost_system -lboost_program_options  $(shell pkg-config --libs-only-l jellyfish-2.0)

all: bin/extract-reads
bin/extract-reads: extract-reads
	cp $< $@

extract-reads: extract-reads.cc
clean:
	rm -f bin/* extract-reads
