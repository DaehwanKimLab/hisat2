# Vendored third-party headers

## sse2neon.h
- Upstream: https://github.com/DLTcollab/sse2neon
- Version: v1.9.1 (pinned)
- License: MIT (see the license header at the top of `sse2neon.h`)
- Purpose: translates the x86 SSE2 intrinsics used by HISAT2's vectorized
  Smith-Waterman implementation (`aligner_swsse_*`) to ARM NEON, enabling
  native aarch64 / Apple Silicon builds. It is placed on the include path
  and included (via `#include "sse2neon.h"`) only on ARM targets; x86 builds
  use the system `<emmintrin.h>` unchanged.

To update: replace the file with a newer pinned release from upstream and bump
the version noted above.

## cpuid.h
- x86-only CPUID helper included by `processor_support.h` on x86 builds to
  detect POPCNT support. Not used on ARM.
