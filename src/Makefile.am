AM_CPPFLAGS = -I$(top_srcdir)/src/ -I$(top_srcdir) $(HTSLIB_CFLAGS)
AM_LDFLAGS = $(HTSLIB_LIBS)

bin_PROGRAMS = chips

chips_SOURCES = main.cpp simulate_reads_main.cpp learn_main.cpp \
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

chips_CPPFLAGS = $(AM_CPPFLAGS) $(AM_PROG_CC_C_O)
chips_CFLAGS = $(CFLAGS)	# Change to AM_CXXFLAGS For -o0 (Valgrind)
chips_CXXFLAGS = $(CXXFLAGS) 	# Change to AM_CXXFLAGS For -o0 (Valgrind)
chips_LDFLAGS = $(AM_LDFLAGS) $(LT_LDFLAGS)
chips_LDADD = $(AM_LDFLAGS) $(LT_LDFLAGS)
