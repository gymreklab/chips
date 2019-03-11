EXTRA_DIST = reconf configure \
	config/git-version-gen \
	$(top_srcdir)/.version

SUBDIRS = m4 src
ACLOCAL_AMFLAGS = -I m4

BUILT_SOURCES = $(top_srcdir)/.version

$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
