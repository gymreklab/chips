AM_CPPFLAGS = -I$(top_srcdir)/src/ -I$(top_srcdir) $(HTSLIB_CFLAGS)
AM_LDFLAGS = $(HTSLIB_LIBS)

bin_PROGRAMS = chipmunk

noinst_LIBRARIES = libchipmunk.a

libchipmunk_a_CPPFLAGS = $(AM_CPPFLAGS) $(HTSLIB_CFLAGS)
libchipmunk_a_SOURCES = simulate_reads_main.cpp learn_main.cpp \
	multithread.h multithread.cpp \
	fragment.h fragment.cpp \
	library_constructor.h library_constructor.cpp \
	options.h options.cpp \
	pulldown.h pulldown.cpp \
	sequencer.h sequencer.cpp \
	common.h common.cpp \
	bam_io.h bam_io.cpp \
	stringops.h stringops.cpp \
	bingenerator.h bingenerator.cpp \
	peak_intervals.h peak_intervals.cpp \
	peak_io_toolbox.h peak_io_toolbox.cpp\
	peak_loader.h peak_loader.cpp\
	ref_genome.h ref_genome.cpp \
	model.h model.cpp json.hpp

chipmunk_SOURCES = main.cpp

chipmunk_LDFLAGS = $(AM_LDFLAGS)
chipmunk_CPPFLAGS = $(AM_CPPFLAGS)
chipmunk_CFLAGS = $(AM_CFLAGS)
chipmunk_CXXFLAGS = $(AM_CSSFLAGS)
chipmunk_LDADD = libchipmunk.a
