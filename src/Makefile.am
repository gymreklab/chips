AM_LDFLAGS = $(HTSLIB_LIBS)

bin_PROGRAMS = asimon

noinst_LIBRARIES = libasimon.a

libasimon_a_CPPFLAGS = $(AM_CPPFLAGS) $(HTSLIB_CFLAGS)
libasimon_a_SOURCES = simulate_reads_main.cpp learn_main.cpp \
	multithread.h multithread.cpp\
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
	ref_genome.h ref_genome.cpp

asimon_SOURCES = main.cpp

asimon_LDFLAGS = $(AM_LDFLAGS)
asimon_CPPFLAGS = $(AM_CPPFLAGS)
asimon_CFLAGS = $(AM_CFLAGS)
asimon_CXXFLAGS = $(AM_CSSFLAGS)
asimon_LDADD = libasimon.a
