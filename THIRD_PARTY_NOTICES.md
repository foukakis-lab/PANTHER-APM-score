# Notices for distributed third-party material

## renv bootstrap

`renv/activate.R` is the renv 1.2.4 bootstrap, retained byte-for-byte from the
validated staging source. Its content matches the installed renv bootstrap
template after the existing version, MD5, and SHA placeholders are substituted.
The generated `.Rprofile` and renv configuration files are also preserved.
No first-party PANTHER license header has been added to these files.

The local renv 1.2.4 DESCRIPTION declares `MIT + file LICENSE`; its LICENSE
records `YEAR: 2026` and `COPYRIGHT HOLDER: Posit Software, PBC`. The following
notice applies to the distributed renv bootstrap, independently of the pending
PANTHER copyright notice.

```text
MIT License

Copyright (c) 2026 Posit Software, PBC

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

## Scope of this notice

No package source trees, installed libraries, model objects, expression matrices,
or participant data are bundled here. R namespace calls and package data loaded
at execution time are external dependencies, not copied package implementations.
Those dependencies retain their own licenses and are not relicensed by this
repository. The inherited lockfile and environment records retain third-party
attribution and provenance metadata. A complete external dependency license
inventory is outside the scope of this distributed-material notice.

The MIT terms above use the standard text from the
[Open Source Initiative](https://opensource.org/license/mit), with the renv
copyright notice taken from the installed package's own LICENSE metadata.
