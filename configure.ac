AC_PREREQ([2.69])
AC_INIT([MELA],[2.0.2],[valerio.bertone@cern.ch],[mela])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AM_MAINTAINER_MODE([disable])

## Set Fortran compiler behaviour
if test "x$FCFLAGS" == "x"; then
  FCFLAGS="-O3"
fi
# Try to respect users' Fortran compiler variables
if test "x$FC" == "x"; then
  if test "x$F77" == "x"; then
    FC="$GFORTRAN"
  else
    FC="$F77"
  fi
fi
AC_PROG_FC
AC_PROG_CXX
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported:\
 use a different --prefix argument])
fi

## Hack to make libtool use the FC Fortran compiler via F77, since
## libtool < 2.x doesn't understand --tag=FC
F77=$FC
AM_FFLAGS=$AM_FCFLAGS
FFLAGS="$FFLAGS $FCFLAGS"
AC_SUBST(FCFLAGS)
AC_SUBST(FFLAGS)
AC_SUBST(F77)
AC_SUBST(AM_FFLAGS)

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

# Check host
AC_CANONICAL_HOST
case $host_os in
  darwin*)
	AM_CPPFLAGS=["$AM_CPPFLAGS -DDARWIN"]
        ;;
esac

## LHAPDF
#AC_ARG_ENABLE([lhapdf],
#  [AC_HELP_STRING(--disable-lhapdf, [don't link to LHAPDF])],
#  [], [enable_lhapdf=yes])
#if test x$enable_lhapdf == xyes; then
#  AC_SEARCH_LHAPDF  
#else
#  AM_CPPFLAGS=["$AM_CPPFLAGS -DNOLHAPDF"]  
#fi
#AM_CONDITIONAL(ENABLE_LHAPDF, [test x$enable_lhapdf == xyes])

## APFEL
#AC_ARG_ENABLE([apfel],
#  [AC_HELP_STRING(--disable-apfel, [don't link to APFEL])],
#  [], [enable_apfel=yes])
#if test x$enable_apfel == xyes; then
#  AC_SEARCH_APFEL  
#else
#  AM_CPPFLAGS=["$AM_CPPFLAGS -DNOAPFEL"]  
#fi
#AM_CONDITIONAL(ENABLE_APFEL, [test x$enable_apfel == xyes])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([lib/Makefile])
AC_CONFIG_FILES([bin/Makefile bin/mela-config])
AC_CONFIG_FILES([run/Makefile])

AC_OUTPUT
