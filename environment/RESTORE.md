# Restoring the recorded environment

Reconstruction and one complete QC/SSPBC validation passed on Windows x64,
R 4.6.0, renv 1.2.4, and Bioconductor 3.23 using the authorized PANTHER inputs.
All 117 previously locked versions were retained; otel 0.2.0 is explicitly
recorded. This establishes the tested Windows/R scope, not other operating
systems or complete reproduction of the paper's APM analyses. The historical
sessionInfo.txt remains the original session record.

Use a fresh, short, non-synced directory outside the repository for the library,
source/binary downloads, cache, sandbox, and build temporary files. Start the
recorded R executable with `--vanilla`. Before startup, set `R_LIBS` and
`R_LIBS_USER` to the new library, `R_LIBS_SITE` to a new empty directory, and
`TMPDIR`, `TMP`, and `TEMP` to the new temporary directory. Set
`RENV_PATHS_ROOT`, `RENV_PATHS_LIBRARY`, `RENV_PATHS_CACHE`, `RENV_PATHS_SOURCE`,
`RENV_PATHS_BINARY`, and `RENV_PATHS_SANDBOX` to separate locations within that
fresh workspace. Do not reuse an installed package library or global cache.

Set these renv environment options before loading renv: `CACHE_ENABLED=FALSE`,
`INSTALL_SHORTCUTS=FALSE`, `USER_LIBRARY=FALSE`, `USER_ENVIRON=FALSE`,
`USER_PROFILE=FALSE`, `AUTO_SNAPSHOT=FALSE`, `CRANDB_ENABLED=FALSE`,
`PAK_ENABLED=FALSE`, `PPM_ENABLED=FALSE`, and `SANDBOX_ENABLED=TRUE`, each
prefixed with `RENV_CONFIG_`. Keep external libraries empty. Disable telemetry
with `OTEL_SDK_DISABLED=true` and set `OTEL_TRACES_EXPORTER`,
`OTEL_METRICS_EXPORTER`, and `OTEL_LOGS_EXPORTER` to `none`.

Bootstrap only renv 1.2.4 into the new library from its exact official CRAN
archive, then load the test project. The recorded source SHA-256 is in
dependency_provenance.tsv. If CRAN has moved that version to its Archive
directory, retrieve that same version and verify the checksum. Do not select
a newer release. Source builds require an already available compatible
toolchain; this validation used the pre-existing Rtools45 GCC 14.3.0 configured
by this R 4.6.0 installation.

With `project` set to the code-only project copy and `lib` to the fresh library:

```r
stopifnot(as.character(getRversion()) == "4.6.0")
.libPaths(c(lib, .Library), include.site = FALSE)
stopifnot(as.character(packageVersion("renv")) == "1.2.4")
renv::load(project = project)
.libPaths(c(lib, renv::paths$sandbox()), include.site = FALSE)
lock <- file.path(project, "renv.lock")
# Complete compiled/binary dependencies before the two URL source packages.
renv::restore(project = project, library = lib, lockfile = lock,
              exclude = c("org.Hs.eg.db", "sspbc"),
              strict = TRUE, retry = FALSE, prompt = FALSE,
              transactional = FALSE)
renv::restore(project = project, library = lib, lockfile = lock,
              packages = c("org.Hs.eg.db", "sspbc"),
              strict = TRUE, retry = FALSE, prompt = FALSE,
              transactional = FALSE)
renv::status(project = project, library = .libPaths(), lockfile = lock)
```

The validation began with an empty library plus fresh renv bootstrap. Initial
URL dependency-order failures were resolved by retaining archive DESCRIPTION
dependencies and aligning 21 Bioconductor repository/remote-type records with
the official sources actually restored. Continuation then installed only the
remaining annotation package. This was not an uninterrupted first-attempt
restore. Do not retry failures using newer versions or install all Suggests.
Base/recommended packages from the recorded R installation may be supplied
through the renv sandbox; codetools 0.2-20 was supplied this way.

The lockfile selects a commit-specific StaafLab sspbc 1.0 archive and the exact
official Bioconductor 3.23 org.Hs.eg.db 3.23.1 archive. Functions/models and the
annotation database matched the validated installation. The original historical
archive identities remain unproven. Selected URLs, hashes, and retrieval dates
are recorded separately from historical installed-package metadata.

Run the unchanged module through a private launcher that selects only the
restored library and sandbox in every R process. Set the documented module
environment variables, keep outputs private, and pass the new Stage 1 output
to Stage 2 and the new final inclusion list to Stage 3. Check every loaded
namespace's version and path; reject original user/site-library packages.

Static dependency scanning does not detect every retained runtime/input
package. renv::status may mark org.Hs.eg.db, otel, tximport, and yaml as recorded
but unused; this is not a version/source mismatch. otel 0.2.0 was loaded during
validated SSPBC execution and is listed under DBI Suggests, not established as
a mandatory DBI Imports dependency. No telemetry exporter is required.

For a future reviewed snapshot, retain the existing explicit package set:

```r
keep <- union(names(renv::lockfile_read(lock)$Packages), "otel")
renv::snapshot(project = project, library = .libPaths(), packages = keep,
               lockfile = "renv.proposed.lock", prompt = FALSE)
```

Review all proposed versions and source records before replacing the lockfile;
retain otel 0.2.0 and the pinned URL sources. Do not snapshot the entire user
library or silently drop packages because static scanning misses them.
See the official renv [record](https://rstudio.github.io/renv/reference/record.html),
[restore](https://rstudio.github.io/renv/reference/restore.html), and
[configuration](https://rstudio.github.io/renv/reference/config.html) documentation.
