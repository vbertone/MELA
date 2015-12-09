# AC_SEARCH_APFEL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_APFEL], [

AC_PATH_PROG(APFELCONFIG, apfel-config, [], [$PATH])
if test -f "$APFELCONFIG"; then
  APFEL_CPPFLAGS=`$APFELCONFIG --cppflags`
  APFEL_LDFLAGS=`$APFELCONFIG --ldflags`
else
  AC_MSG_ERROR([APFEL cannot be found!])
  exit 1
fi
AC_SUBST(APFEL_CPPFLAGS)
AC_SUBST(APFEL_LDFLAGS)
$1
])
