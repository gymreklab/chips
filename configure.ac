AC_INIT([asimon], [1.0], [mgymrek@ucsd.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])
AC_PROG_RANLIB
AM_PROG_AR
AC_PROG_CXX

# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

# Custom Flags to add to CFLAGS/CXXFLAGS
asimon_CFLAGS="-D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\""
CFLAGS="$CFLAGS $asimon_CFLAGS"
CXXFLAGS="$CXXFLAGS $asimon_CFLAGS"

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
