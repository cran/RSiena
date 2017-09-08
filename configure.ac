# vim:fdm=marker:

# Based on the autoconf.ac from the Rmpi package
# <http://www.stats.uwo.ca/faculty/yu/Rmpi/>. Run `autoconf` to produce the
# 'configure' script.

# - Run with automatic detection: `R CMD INSTALL {pkg}`

# - Run with a prefix folder. Subfolders will be searched for binarys.
#   `R CMD INSTALL {pkg} --configure-args=--with-mpi=/path/to/location`

# - Fully specify paths and type.
#   `R CMD INSTALL {pkg} --configure-args="        \\
#    --with-Rmpi-include=/path/to/mpi_include_dir  \\
#    --with-Rmpi-libpath=/path/to/mpi_lib_dir      \\
#    --with-Rmpi-type={OPENMPI,MPICH,MPICH2,LAM}"`

AC_PREREQ([2.69])
AC_INIT([RSienaTest], 1.1-289)
AC_CONFIG_SRCDIR([src])

# Checks (autoscan)                                                       {{{1
AC_PROG_CXX
AC_PROG_CC
# Checks for header files
AC_CHECK_HEADERS([sys/time.h unistd.h])
# Checks for typedefs, structures, and compiler characteristics
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])
# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime gettimeofday memmove memset pow select sqrt])

# Args (mpi discovery)                                                    {{{1
dnl --with-mpi-include-path
MPI=yes
AC_ARG_WITH([mpi-include-path],
            AS_HELP_STRING([--with-mpi-include-path=DIR], [location of MPI header files]),
            [if test "${withval}" != "no" ; then
             MPI_INCLUDE_PATH="${withval}"
             if test ! -f "${MPI_INCLUDE_PATH}/mpi.h"; then
               AC_MSG_ERROR([Value of --with-include-path does not contain mpi.h])
               fi
             else
               if test -n "${MPI_INCLUDE_PATH}" ; then
                 MPI_INCLUDE_PATH="${MPI_INCLUDE_PATH}"
                 if test ! -f "${MPI_INCLUDE_PATH}/mpi.h"; then
                   AC_MSG_ERROR([Value of MPI_INCLUDE_PATH does not contain mpi.h])
                   fi
                   echo "Setting MPI include path from MPI_INCLUDE_PATH"
                 fi
               fi])

dnl --with-mpi-lib-path
AC_ARG_WITH([mpi-lib-path],
            AS_HELP_STRING([--with-mpi-lib-path=DIR], [location of MPI library files]),
            [if test "${withval}" != "no" ; then
             MPI_LIB_PATH="${withval}"
           else
             if test -n "${MPI_LIB_PATH}" ; then
               MPI_LIB_PATH="${MPI_LIB_PATH}"
               echo "Setting MPI lib path from MPI_LIB_PATH"
             fi
           fi])

if test -n "${MPI_INCLUDE_PATH}" -a -n "${MPI_LIB_PATH}" ; then
  MPI_INCLUDES="-I${MPI_INCLUDE_PATH}"
  MPI_LIBS="-L${MPI_LIB_PATH}"
elif test -n "${MPI_INCLUDE_PATH}" -o -n "${MPI_LIB_PATH}" ; then
  AC_MSG_ERROR([Must specify both the include and lib path])
fi

dnl --with-mpi-type
AC_ARG_WITH([mpi-type],
            AS_HELP_STRING([--with-mpi-type=MPI_TYPE], [the type of MPI: OPENMPI, LAM or MPICH]),
            [if test "${withval}" != "no" ; then
              MPI_TYPE="${withval}"
            else
              if test -n "${MPI_TYPE}" ; then
                MPI_TYPE="${MPI_TYPE}"
                echo "Setting MPI type from MPI_TYPE"
              fi
            fi])

dnl begin auto configure paths
if test -z "${MPI_INCLUDE_PATH}" ; then

  dnl --with-mpi-root
  AC_ARG_WITH(mpi-root,
              AS_HELP_STRING([--with-mpi-root=DIR], [location of top-level MPI directory]),
              [if test "${withval}" != "no" ; then
                MPI_ROOT="${withval}"
                MPI_INCLUDE_PATH="${MPI_ROOT}/include"
                MPI_LIB_PATH="${MPI_ROOT}/lib"
              fi])

  if test -z "$MPI_ROOT" ; then
    dnl auto configure MPI_ROOT, (MPI_TYPE, MPI_INCLUDE_PATH)
    AC_MSG_CHECKING([for mpi root directory])
    for d in /opt/lib /usr/lib /usr /usr/local/lib /usr/local \
      /usr/lib64/mpi/gcc/openmpi; do
      if test -f $d/include/mpi.h && test -d $d/lib/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        break
      elif test -f $d/include/openmpi/mpi.h && test -d $d/lib/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/include/openmpi/mpi.h && test -d $d/lib64/openmpi; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/openmpi/include/mpi.h && test -d $d/openmpi/lib; then
        MPI_ROOT=$d/openmpi
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include/openmpi
        break
      elif test -f $d/include/mpi.h && test -d $d/lib64; then
        MPI_ROOT=$d
        MPI_TYPE="OPENMPI"
        MPI_INCLUDE_PATH=$d/include
        break
      elif test -f $d/lam/include/mpi.h && test -f $d/lam/lib/libmpi.so; then
        MPI_ROOT=$d/lam
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/include/lam/mpi.h && test -f $d/lib/lam/libmpi.so; then
        MPI_ROOT=$d
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/include/lam/mpi.h && test -f $d/lib64/lam/libmpi.so; then
        MPI_ROOT=$d
        MPI_TYPE="LAM"
        MPI_INCLUDE_PATH=$d/include/lam
        break
      elif test -f $d/mpich/include/mpi.h; then
        MPI_ROOT=$d/mpich
        MPI_TYPE="MPICH"
        MPI_INCLUDE_PATH=$d/mpich/include
        break
      elif test -f $d/mpi/include/mpi.h; then
        MPI_ROOT=$d/mpi
        MPI_INCLUDE_PATH=$d/mpi/include
        break
      elif test -f $d/include/mpi.h; then
        MPI_ROOT=$d
        break
      fi
    done
    AC_MSG_RESULT([$MPI_ROOT])
  fi

  if test -n "$MPI_ROOT"; then
    dnl auto configure MPI_INCLUDE_PATH
    AC_MSG_CHECKING([for mpi.h])
    for d in include include/openmpi include/lam ; do
      if test -f ${MPI_ROOT}/$d/mpi.h ; then
        AC_MSG_RESULT([${MPI_ROOT}/$d])
        MPI_INCLUDE_PATH="${MPI_ROOT}/$d"
        break
      fi
    done
    if test -z "$MPI_INCLUDE_PATH" ; then
      AC_MSG_RESULT([no])
      MPI=no
    fi
  else
    AC_CHECK_HEADER(mpi.h, , [MPI=no])
  fi
  if test -n "$MPI_INCLUDE_PATH" ; then
    MPI_INCLUDES="-I${MPI_INCLUDE_PATH}"
  fi

  dnl auto configure MPI_LIB_PATH, (MPI_INCLUDES)
  AC_MSG_CHECKING([for libmpi.{a,so,dylib}])
  if test -f ${MPI_ROOT}/lib/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  elif test -f ${MPI_ROOT}/libmpi.so; then
    dnl Fedora 10
    MPI_LIB_PATH="${MPI_ROOT}"
  elif test -f ${MPI_ROOT}/lib/openmpi/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib/openmpi"
  elif test -f ${MPI_ROOT}/lib64/openmpi/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64/openmpi"
  elif test -f ${MPI_ROOT}/lib64/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64"
  elif test -f ${MPI_ROOT}/lib/lam/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib/lam"
    MPI_INCLUDES="$MPI_INCLUDES -I${MPI_INCLUDE_PATH}/32"
  elif test -f ${MPI_ROOT}/lib64/lam/libmpi.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64/lam"
    MPI_INCLUDES="$MPI_INCLUDES -I${MPI_INCLUDE_PATH}/64"
  elif test -f ${MPI_ROOT}/lib/libmpich.so; then
    MPI_LIB_PATH="${MPI_ROOT}/lib"
    AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],,[#include <mpi.h>])
  elif test -f ${MPI_ROOT}/lib64/libmpich.a; then
    MPI_LIB_PATH="${MPI_ROOT}/lib64"
    AC_CHECK_DECLS([MPICH2, MPICH2],[  MPI_DEFS="-DMPI2" ],,[#include <mpi.h>])
  elif test -f ${MPI_ROOT}/lib/libmpi.dylib; then
    dnl OSX (10.6)
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  else
    AC_CHECK_LIB(mpi, main, [MPI_LIB_PATH="${MPI_ROOT}/lib"], )
  fi
  if test -n "$MPI_LIB_PATH" ; then
    AC_MSG_RESULT([$MPI_LIB_PATH])
  fi
fi
dnl end auto configure paths
if test -z $MPI_TYPE ; then
  if   test -d ${MPI_ROOT}/openmpi;          then MPI_TYPE="OPENMPI"
  elif test -d ${MPI_ROOT}/lib/openmpi;      then MPI_TYPE="OPENMPI"
  elif test -d ${MPI_ROOT}/lib64/openmpi;    then MPI_TYPE="OPENMPI"
  elif test -f ${MPI_ROOT}/lib/liblam.so;    then MPI_TYPE="LAM"
  elif test -d ${MPI_ROOT}/lib/lam;          then MPI_TYPE="LAM"
  elif test -d ${MPI_ROOT}/lib64/lam;        then MPI_TYPE="LAM"
  elif test -f ${MPI_ROOT}/lib/libmpich.a;   then MPI_TYPE="MPICH"
  elif test -f ${MPI_ROOT}/lib64/libmpich.a; then MPI_TYPE="MPICH"
  fi
fi

if test "$MPI_TYPE" = "OPENMPI" ; then
  AC_CHECK_PROG(ORTED, orted, yes, no)
  if test "$ORTED" = no ; then
    AC_MSG_NOTICE([Cannot find orted. Rmpi needs orted to run.])
  fi
fi

if test "$MPI_TYPE" = "LAM" ; then
  echo "Try to find liblam.so ..."
  if test -f ${MPI_ROOT}/lib/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib"
    MPI_LIB_PATH="${MPI_ROOT}/lib"
  elif test -f ${MPI_ROOT}/lib/lam/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib/lam"
    MPI_LIB_PATH="${MPI_ROOT}/lib/lam"
  elif test -f ${MPI_ROOT}/lib64/lam/liblam.so; then
    echo "Found liblam in ${MPI_ROOT}/lib64/lam"
    MPI_LIB_PATH="${MPI_ROOT}/lib64/lam"
  else
    AC_CHECK_LIB(lam, main, [ echo "liblam not found. Probably not LAM-MPI" ])
  fi
fi

if test -n "$MPI_LIB_PATH" ; then
  MPI_LIBS="-L${MPI_LIB_PATH}"
else
  MPI=no
fi

AC_CHECK_LIB(util, openpty, [ MPI_LIBS="$MPI_LIBS -lutil" ])
AC_CHECK_LIB(pthread, main, [ MPI_LIBS="$MPI_LIBS -lpthread" ])

MPI_DEFS="-D$MPI_TYPE"
if test "$MPI_TYPE" = "OPENMPI" -o "$MPI_TYPE" = "LAM" -o "$MPI_TYPE" = "MPICH2"; then
  MPI_DEFS="$MPI_DEFS -DMPI2"
fi
OSTYPE="`uname`"
if test "$OSTYPE" = "Darwin"; then
  MPI_DEFS="$MPI_DEFS -DMAC"
fi

PKG_CPPFLAGS="$MPI_INCLUDES $MPI_DEFS"
if test "$MPI" = "yes" ; then
  case "$MPI_TYPE" in
    OPENMPI) PKG_LIBS="-lmpi_cxx -lmpi ${MPI_LIBS}";;
    LAM)     PKG_LIBS="-lmpi_cxx -llam ${MPI_LIBS}";;
    MPICH)   PKG_LIBS="-lmpich -lmpl -lopa ${MPI_LIBS}";;
    MPICH2)  PKG_LIBS="-lmpich -lmpl -lopa ${MPI_LIBS}";;
    *)
      AC_MSG_NOTICE([Unknown type of mpi, compiling without mpi support. Use --with-mpi-xxxx to specify it.])
      PKG_CPPFLAGS=
      PKG_LIBS=
      ;;
  esac
else
  PKG_CPPFLAGS=
  PKG_LIBS=
fi

PKG_SOURCES=$(cat src/sources.list)

AC_SUBST(PKG_SOURCES)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)

AC_OUTPUT(src/Makevars)
