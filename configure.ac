AC_INIT([tulip],
	[m4_esyscmd([config/git-version-gen .tarball-version])], 
	[Melissa Gymrek mgymrek@eng.ucsd.edu],
	[tulip])
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
PKG_CHECK_MODULES([HTSLIB],[htslib])


AC_SUBST(LT_LDFLAGS)

# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

# Compiler warnings to add
#COMPILER_WARNINGS="-Wextra -Wswitch-default -Wno-strict-aliasing"

# Check for PTHREADS
AX_PTHREAD([],[AC_MSG_ERROR([Pthread support required for tulip])])

# Custom Flags to add to CFLAGS/CXXFLAGS
TULIP_CFLAGS="-D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\" -std=gnu++11"

# Set the final value for CFLAGS/CXXFLAGS
CFLAGS="$CFLAGS $COMPILER_WARNINGS $PTHREAD_CFLAGS $TULIP_CFLAGS"
CXXFLAGS="$CXXFLAGS $COMPILER_WARNINGS $PTHREAD_CFLAGS $TULIP_CFLAGS"

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
AC_MSG_RESULT([    PTHREAD_LIBS:  $PTHREAD_LIBS])
AC_MSG_RESULT([    PTHREAD_CFLAGS: $PTHREAD_CFLAGS])
AC_MSG_RESULT([    TULIP_CFLAGS: $TULIP_CFLAGS])
AC_MSG_RESULT([])
AC_MSG_RESULT([ Default installtion directories:])
AC_MSG_RESULT([    programs: ${prefix}/bin/ ])
AC_MSG_RESULT([])
AC_MSG_RESULT([ To change installtion path, use:])
AC_MSG_RESULT([    ./configure --prefix NEW-PATH])
AC_MSG_RESULT([])
