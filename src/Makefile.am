bin_PROGRAMS = asimon

noinst_LIBRARIES = libasimon.a

libasimon_a_CPPFLAGS = $(AM_CPPFLAGS)
libasimon_a_SOURCES = simulate_reads_main.cpp learn_main.cpp \
	fragment.h fragment.cpp \
	library_constructor.h library_constructor.cpp \
	options.h options.cpp \
	pulldown.h pulldown.cpp \
	sequencer.h sequencer.cpp \
	shearer.h shearer.cpp

asimon_SOURCES = main.cpp

asimon_LDFLAGS = $(AM_LDFLAGS)
asimon_CPPFLAGS = $(AM_CPPFLAGS)
asimon_CFLAGS = $(AM_CFLAGS)
asimon_CXXFLAGS = $(AM_CSSFLAGS)
asimon_LDADD = libasimon.a
