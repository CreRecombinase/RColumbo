AC_INIT([RColumbo], 0.0.1, [nwknoblauch@gmail.com])

dnl Now find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

AC_CHECK_PROG(H5CPP,h5c++,`which h5c++`)
if test -z "${H5CPP}"; then
  AC_MSG_ERROR(['h5c++' does not seem to be installed on your platform.
  	Please install the hdf5 library with C++ interface enabled (--enable-cxx=yes).
  	The required HDF5 library files can be installed as follows:
  	    - Debian-based (e.g. Debian >= 8.0, Ubuntu >= 15.04): 'sudo apt-get install libhdf5-dev'
        - Old Debian-based (e.g Debian < 8.0, Ubuntu < 15.04): Install from source (see INSTALL)
        - OS X using Homebrew: 'brew install homebrew/science/hdf5 --enable-cxx'
        - RPM-based (e.g Fedora): 'sudo yum install hdf5-devel'])
fi

if test -z "${LIBS}"; then
  dnl prefer ggrep if installed
  AC_CHECK_PROG(GREP,ggrep,`which ggrep`, `which grep`)
  
  dnl check if $GREP supports option -o
  AC_MSG_CHECKING([whether ${GREP} accepts -o])
  if `echo '' | ${GREP} -o '' >/dev/null 2>&1`
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([${GREP} does not support option -o])
  fi
  [
  LIBS=`${H5CPP} -show | ${GREP} -o "\-L[^[:space:]]\+" | tr "\n" " "` 
  ]
  
  dnl Fix for newer Debian versions
  LIBS="${LIBS} -L. -L${HOME}/lib -L${HOME}/lib64 -lhdf5_cpp -lhdf5 -lz -lblosc"

fi

if test -z "${CPPFLAGS}"; then
  dnl prefer ggrep if installed
  AC_CHECK_PROG(GREP,ggrep,`which ggrep`, `which grep`)

  dnl check if $GREP supports option -o
  AC_MSG_CHECKING([whether ${GREP} accepts -o])
  if `echo '' | ${GREP} -o '' >/dev/null 2>&1`
  then
    AC_MSG_RESULT([yes])
  else
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([${GREP} does not support option -o])
  fi

  echo 'int main() { return; }' > myconftest.cpp

  [
  CPPFLAGS=`${H5CPP} -show -c myconftest.cpp | ${GREP} -o "\-I[^[:space:]]\+" | tr "\n" " "`
  CPPFLAGS="${CPPFLAGS} -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -DARMA_64BIT_WORD -I${HOME}/include"
  ]

  rm -f myconftest.cpp
fi


# Checks for header files.
dnl AC_CHECK_HEADERS([inttypes.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
dnl AC_CHECK_HEADER_STDBOOL
dnl AC_C_INLINE
dnl AC_TYPE_INT16_T
dnl AC_TYPE_INT32_T
dnl AC_TYPE_INT64_T
dnl AC_TYPE_INT8_T
dnl AC_TYPE_SIZE_T
dnl AC_TYPE_UINT16_T
dnl AC_TYPE_UINT32_T
dnl AC_TYPE_UINT64_T
dnl AC_TYPE_UINT8_T

# Checks for library functions.
dnl AC_FUNC_MALLOC
dnl AC_CHECK_FUNCS([floor memset sqrt strdup strerror strtol])


dnl substitute CPPFLAGS and LIBS
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
dnl and do subsitution in the src/Makevars.in 
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
