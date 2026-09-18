# AI Disclosure

AI tools (Claude, Anthropic) were used throughout this project as a force
multiplier for a human-directed research workflow — not as a substitute for
the researcher's own design, judgment, or scientific reasoning. All code
that AI tools touched was built on top of human-written analysis code and
methodology; AI did not originate the science.

Specific uses include:

- **Organization** — structuring repositories, file layout, and documentation
- **Code streamlining** — refactoring and cleaning up human-written code for
  clarity and reuse, without changing its underlying logic
- **Version control** — commit hygiene and git workflow management
- **Environment building** — consolidating and managing software dependencies
  (e.g., conda/renv environments, pipeline configuration)
- **Troubleshooting** — debugging pipeline errors, dependency conflicts, and
  unexpected tool behavior
- **Designing pilot runs** — scaffolding small-scale test/validation runs
  ahead of full analyses

Experimental design, data interpretation, and all scientific conclusions in
this repository are the author's own.

**`nextflow_pipeline/` specifically**: its Nextflow structure (modules,
subworkflows, parameterization, conda environments, test scaffolding) was
written by Claude and reviewed/validated by the author. The bulk RNA-seq and
scRNA-seq/velocity logic it runs is a direct, parameterized port of this
repo's own author-written `bulkSeq_code.R`/`seurat_code.R`/`scvelo_code.py`.
Its `scenic` mode ports a separate pipeline (originally Claude-written,
author-reviewed) from an unrelated project and does not apply to this
repo's data — see `nextflow_pipeline/README.md`.
