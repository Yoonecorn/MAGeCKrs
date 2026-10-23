# Third-Party Notices

This project uses and adapts data and components from the original [MAGeCK2](https://bitbucket.org/liulab/mageck2) software.

## MAGeCK2

The embedded Normal Quantile table (`src/mageck_test/qnorm_table.txt`) and some numeric validations are derived directly from MAGeCK2, specifically from reference commit `3616dd11e10fcd8d586d6a8eb76b2faa493d6aae`, which is licensed under the BSD 3-Clause License.

```
BSD 3-Clause License

Copyright (c) 2021, Wei Li
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Cephes normal-tail implementation

The deterministic Rust translation of Cephes Math Library Release 2.0 `ndtr.c`
in `src/mle/wald.rs` preserves this source notice:

```text
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier.

Permission to use, copy, modify, and distribute this software is freely
granted, provided that this notice is preserved.
```

## Locked nalgebra published-source records

The following exact published files are preserved under `docs/third-party/`.
Their `.crate` archive SHA-256 values match the corresponding `Cargo.lock`
checksums; this authenticates the observed records to the locked package bytes.
It does not resolve conflicting terms or select one license as authoritative.
The project MIT license does not replace or relicense these third-party records.

| Locked package | Published manifest declaration | Included LICENSE | Narrow unresolved issue |
|---|---|---|---|
| nalgebra 0.32.6 | `BSD-3-Clause` in both manifest forms | Apache License 2.0 text with Copyright 2020 Sébastien Crozet | The published declaration and included text disagree; no choice or legal reconciliation is asserted. |
| nalgebra-macros 0.2.2 | `Apache-2.0` in both manifest forms | Literal `../LICENSE`, 10 bytes, no final newline; a regular file | The published archive does not resolve this workspace pointer to a full license file. |

Preserved files: [nalgebra Cargo.toml](docs/third-party/nalgebra-0.32.6/Cargo.toml),
[Cargo.toml.orig](docs/third-party/nalgebra-0.32.6/Cargo.toml.orig),
[LICENSE](docs/third-party/nalgebra-0.32.6/LICENSE), and
[nalgebra-macros Cargo.toml](docs/third-party/nalgebra-macros-0.2.2/Cargo.toml),
[Cargo.toml.orig](docs/third-party/nalgebra-macros-0.2.2/Cargo.toml.orig),
[LICENSE](docs/third-party/nalgebra-macros-0.2.2/LICENSE).
The nalgebra LICENSE retains its published CRLF bytes. The nalgebra-macros
pointer remains exact; it is not a symlink and no generic text is silently
substituted under that filename. Published third-party manifest URLs are retained
as source records; they are not a MAGeCKrs repository URL.

| Archive | Bytes | SHA-256 (also Cargo.lock checksum) |
|---|---:|---|
| nalgebra-0.32.6.crate | 375297 | `7b5c17de023a86f59ed79891b2e5d5a94c705dbe904a5b5c9c952ea6221b03e4` |
| nalgebra-macros-0.2.2.crate | 7545 | `254a5372af8fc138e36684761d3c0cdb758a4410e938babcff1c860ce14ddbfc` |

Authoritative upstream clarification or a reviewed licensing decision is still
needed before claiming those issues are legally resolved. No dependency upgrade
or substitution is performed here. The records above preserve the MAGeCK2 and Cephes notices and these two locked
published-source records; additional dependency disclosures follow below.
This file does not claim complete legal clearance of every dependency.

## Source-only dependency scope

MAGeCKrs distributes its own source tree, without vendoring full dependency
source trees or distributing a compiled executable. `Cargo.lock` records dependency
versions and checksums; Cargo obtains dependencies during a build. This scope does
not remove all notice duties. The MAGeCKrs MIT license does not relicense
third-party packages, and this file does not claim complete external legal clearance.

## Optional Metal dependency disclosures

The publication decision accepts the following evidence-limited disclosures.
Both packages declare MIT and occur in the optional Metal dependency path.
Complete historical notice evidence was not found for these exact releases;
no redistribution prohibition was identified in the reviewed evidence. This is a
publication decision, not legal advice or a claim of complete legal clearance.

- [block 0.1.6](docs/third-party/block-0.1.6/NOTICE.md): the exact archive lacks
  complete release-specific notice text; the reviewed upstream notice work was
  not merged for this release.
- [malloc_buf 0.0.6](docs/third-party/malloc_buf-0.0.6/NOTICE.md): the exact archive
  lacks complete release-specific notice text; a later upstream MIT notice exists,
  but is not treated as proof that it applied to version 0.0.6.

MAGeCKrs does not relicense either package. No external legal clearance is claimed.
