AC_INIT([asimon],
	[m4_esyscmd([config/git-version-gen .tarball-version])], 
	[Melissa Gymrek mgymrek@eng.ucsd.edu],
	[asimon])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/main.cpp])
AM_INIT_AUTOMAKE([foreign -Werror subdir-objects])
AM_PROG_AR

LT_INIT([static])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AM_PROG_CC_C_O 

# Check for required packages
#PKG_CHECK_MODULES([GSL],[gsl])
PKG_CHECK_MODULES([HTSLIB],[htslib])
#PKG_CHECK_MODULES([NLOPT],[nlopt])

# To compile a static executable (before binary packaging?),
# use:
#    ./configure --enable-all-static
#
# NOTE:
#  1.  Using "./configure LDFLAGS=-static" will not work, because
#      libtool will override the LDFLAG files at times.
#  2.  The flag 'LT_LDFLAGS' is used in ./src/Makefile.am with LDADD.
#
# 'all-static' macro copied from subversion's configure.ac
AC_ARG_ENABLE(all-static,
		AS_HELP_STRING([--enable-all-static],
                             [Build completely static (standalone) binaries.]),
                             [
                              if test "$enableval" = "yes" ; then
                                      LT_LDFLAGS="-all-static  $LT_LDFLAGS"
                              elif test "$enableval" != "no" ; then
                                      AC_MSG_ERROR([--enable-all-static doesn't accept argument])
                                      fi
                                      ])

AC_SUBST(LT_LDFLAGS)

# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

# Compiler warnings to add
#COMPILER_WARNINGS="-Wextra -Wswitch-default -Wno-strict-aliasing"

# Custom Flags to add to CFLAGS/CXXFLAGS
ASIMON_CFLAGS="-D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\""

# Set the final value for CFLAGS/CXXFLAGS
CFLAGS="$CFLAGS $COMPILER_WARNINGS $ASIMON_CFLAGS"
CXXFLAGS="$CXXFLAGS $COMPILER_WARNINGS $ASIMON_CFLAGS"

AC_CONFIG_FILES([
   m4/Makefile
   Makefile
   src/Makefile
])

AC_OUTPUT
AC_MSG_RESULT([])
AC_MSG_RESULT([ Configuration summary for $PACKAGE])
AC_MSG_RESULT([    version $VERSION])
AC_MSG_RESULT([    CFLAGS:   $CFLAGS])
AC_MSG_RESULT([    CPPFLAGS: $CPPFLAGS])
AC_MSG_RESULT([    CXXFLAGS: $CXXFLAGS])
AC_MSG_RESULT([    LDFLAGS:  $LDFLAGS])
AC_MSG_RESULT([    HTSLIB_LIBS:  $HTSLIB_LIBS])
AC_MSG_RESULT([    HTSLIB_CFLAGS: $HTSLIB_CFLAGS])
AC_MSG_RESULT([    ASIMON_CFLAGS: $ASIMON_CFLAGS])
AC_MSG_RESULT([])
AC_MSG_RESULT([ Default installtion directories:])
AC_MSG_RESULT([    programs: ${prefix}/bin/ ])
AC_MSG_RESULT([])
AC_MSG_RESULT([ To change installtion path, use:])
AC_MSG_RESULT([    ./configure --prefix NEW-PATH])
AC_MSG_RESULT([])
