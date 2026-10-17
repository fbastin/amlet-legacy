#!/usr/bin/env bash
#
# Build ORATIO -> OPHELIA -> AMLET into a local prefix, without root access.
#
#   ./build.sh                 # installs into ./local
#   PREFIX=$HOME/opt ./build.sh
#   BLAS_LIBS="-lopenblas" ./build.sh      # force a BLAS/LAPACK choice
#
# If the source tree sits on a filesystem mounted "noexec" (e.g. an external
# backup disk), put the build and install directories elsewhere:
#   PREFIX=$HOME/amlet BUILD=/tmp/amlet-build bash build.sh
#
# Requirements: a C compiler (gcc or clang), make, and a BLAS library that
# exports the CBLAS interface (cblas_dgemv, ...) plus LAPACK (dpotrf_).
# On Debian/Ubuntu: apt install build-essential libblas-dev liblapack-dev
#   (or libopenblas-dev, which provides both).
# Bison/flex, autotools and LaTeX are NOT needed: the generated parser,
# scanner and configure scripts are shipped in the source trees.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${PREFIX:-$ROOT/local}"
BUILD="${BUILD:-$ROOT/build}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 2)}"
CC="${CC:-gcc}"
CFLAGS="${CFLAGS:--O2}"

ORATIO=oratio-0.7.5
OPHELIA=ophelia-0.7.4
AMLET=amlet-0.13.1

log() { printf '\n==> %s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

mkdir -p "$PREFIX/lib" "$PREFIX/include" "$BUILD"
PREFIX="$(cd "$PREFIX" && pwd)"
BUILD="$(cd "$BUILD" && pwd)"

for d in "$BUILD" "$PREFIX"; do
  printf '#!/bin/sh\nexit 0\n' > "$d/.exec-test" && chmod +x "$d/.exec-test"
  if ! "$d/.exec-test" 2>/dev/null; then
    rm -f "$d/.exec-test"
    die "$d is on a noexec filesystem; set PREFIX and BUILD to another location (see header of build.sh)"
  fi
  rm -f "$d/.exec-test"
done

# ---------------------------------------------------------------------------
# BLAS / LAPACK detection
# ---------------------------------------------------------------------------
# The historical configure scripts only accept -lopenblas, -lgoto2, -lcblas,
# -lgslcblas+-lblas or -latlas. Modern distributions often ship the CBLAS
# symbols inside libblas.so instead; in that case we expose it to configure
# under the name libcblas.so through a link-time-only shim directory (the
# binaries end up depending on the real library's soname, e.g. libblas.so.3).

can_link() {  # can_link <symbol> <libs...>
  local sym="$1"; shift
  printf 'char %s(void);\nint main(void){return (int)%s();}\n' "$sym" "$sym" \
    | "$CC" -x c - -o "$BUILD/conftest" "$@" >/dev/null 2>&1
}

SHIM="$BUILD/blas-shim"
mkdir -p "$SHIM"

if [ -n "${BLAS_LIBS:-}" ]; then
  can_link cblas_dgemv $BLAS_LIBS || die "BLAS_LIBS='$BLAS_LIBS' does not provide cblas_dgemv"
else
  BLAS_LIBS=""
  for cand in openblas cblas blas; do
    if can_link cblas_dgemv "-l$cand"; then BLAS_LIBS="-l$cand"; break; fi
  done
  [ -n "$BLAS_LIBS" ] || die "no BLAS library with the CBLAS interface found (install libopenblas-dev or libblas-dev)"
fi

case "$BLAS_LIBS" in
  -lopenblas*|-lcblas*) ;;
  *)
    real="$("$CC" -print-file-name="lib${BLAS_LIBS#-l}.so")"
    [ -e "$real" ] || die "cannot locate shared library for $BLAS_LIBS"
    ln -sf "$real" "$SHIM/libcblas.so"
    ;;
esac

if can_link dpotrf_ $BLAS_LIBS; then
  LAPACK_LIBS=""
elif can_link dpotrf_ -llapack $BLAS_LIBS; then
  LAPACK_LIBS="-llapack"
else
  die "no LAPACK found (dpotrf_); install liblapack-dev or libopenblas-dev"
fi
rm -f "$BUILD/conftest"

log "Using BLAS: $BLAS_LIBS   LAPACK: ${LAPACK_LIBS:-(inside BLAS)}"

export CC
export CPPFLAGS="-I$PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$PREFIX/lib -L$SHIM ${LDFLAGS:-}"
export LIBS="$LAPACK_LIBS $BLAS_LIBS ${LIBS:-}"

# git does not keep file dates: after a clone, configure.ac or a .y file may
# look newer than the files generated from it, and make would then try to
# rerun autotools or bison. Make each generated file newer than its sources.
refresh() {  # refresh <generated file> <sources...>
  local out="$1" src; shift
  [ -e "$out" ] || return 0
  for src in "$@"; do
    if [ -e "$src" ] && [ "$src" -nt "$out" ]; then touch "$out"; return 0; fi
  done
}

fix_timestamps() {  # fix_timestamps <tree>
  local d="$1" f
  refresh "$d/aclocal.m4" "$d/configure.ac" "$d"/m4/*.m4 "$d/acinclude.m4"
  refresh "$d/configure" "$d/configure.ac" "$d/aclocal.m4"
  refresh "$d/config.h.in" "$d/configure.ac" "$d/aclocal.m4"
  find "$d" -name Makefile.in | while read -r f; do
    refresh "$f" "${f%.in}.am" "$d/configure.ac" "$d/aclocal.m4"
  done
  find "$d" \( -name '*.y' -o -name '*.l' \) | while read -r f; do
    refresh "${f%.*}.c" "$f"
    refresh "${f%.*}.h" "${f%.*}.c"
  done
}

# configure_build <source dir> <subdirs to build and install...>
# The old Makefiles do not support out-of-tree (VPATH) builds, so each source
# tree is copied into $BUILD and built there; src/ stays untouched. The copy
# is refreshed on every run, so edits in src/ are picked up incrementally
# (remove $BUILD for a clean rebuild).
configure_build() {
  local src="$ROOT/src/$1"; shift
  local dir="$BUILD/$(basename "$src")"
  mkdir -p "$dir"
  cp -pR "$src/." "$dir/"   # -p keeps timestamps; otherwise make reruns automake
  fix_timestamps "$dir"
  log "Configuring $(basename "$src")"
  (cd "$dir" && "${CONFIG_SHELL:-/bin/sh}" ./configure --prefix="$PREFIX" CFLAGS="$CFLAGS" \
     > configure.log 2>&1) || die "configure failed, see $dir/configure.log"
  for sub in "$@"; do
    log "Building $(basename "$src")/$sub"
    (cd "$dir/$sub" && make -j"$JOBS" > make.log 2>&1 && make install >> make.log 2>&1) \
      || die "build failed, see $dir/$sub/make.log"
  done
}

configure_build "$ORATIO"  oratio ctools dstruct ntools stat graph ran sim moca
configure_build "$OPHELIA" ophelia lapack
configure_build "$AMLET"   amlet src

log "Done. Executable: $PREFIX/bin/amlet"
(cd "$BUILD" && "$PREFIX/bin/amlet" -V | head -1) || true   # amlet writes amlet.log in the cwd
