# Changelog

A release that changes any output bit of an existing function names that
function under the release, even for a one-ulp change. `tests/golden.txt`
pins a hash of every function's outputs on four types, and the tests fail when
one moves; `GOLDEN_UPDATE=1 cargo test golden` re-pins it and stamps each
changed row with the crate version, whose section here must then name the
function in backticks.

Releases before 3.1.1 are described on the
[releases page](https://github.com/GeEom/fixed_analytics/releases).

## 3.1.1

No output changes. The package leaves out the toolchain and CI files.
