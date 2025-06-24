
# Linear epitope prediction using AlphaFold3

## Setup

Create a .env file in the project root:

```bash
CONDA_PREFIX=/path/to/conda/installation
```

Swap the path in `envs/env.yaml` to the project root to your cloned path:

```yaml
...
  - pip:
    # swap me for git URL
    - /path/to/cloned/repo
...
```

Install the nextflow and project environments:

```bash
# project environment
conda env create --name <name> --file envs/env.yaml
# nf-core/nextflow env 
conda env create --name <name> --file envs/nf-core.yaml
```

If running on Gemini, ensure that Alphafold MSAs will use your scratch as a tmp directory (some MSA intermediate files are larger than `/tmp` on compute nodes)

```bash
vim ~/.nextflow/config
```

```bash
params.msa_tmpdir = "/scratch/<username>/tmp"
```