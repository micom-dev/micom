# Running MICOM with Snakemake

Several steps of MICOM can be quite time consuming. Luckily, they are also parallelizable! Here we describe how to use the included Snakemake workflow to distribute those tasks to the cloud or HPC, and hand the results back to the tutorial.

## Prerequisites

- an environment containing MICOM and all its dependencies, as well as Snakemake
  -- eg, `mamba create -n micom snakemake python -y ; mamba activate micom ;  pip install numpy==2.1.1 Cython biom-format==2.1.16 micom seaborn`
- a configured Snakemake Profile

## Input Data

This workflow requires three files:

- a long-format sample counts table with columns as follows:
  - `sample_id`
  - `genus`
  - `abundance`
- a MICOM-formatted AGORA model bundle: `wget -O agora103_genus.qza https://zenodo.org/record/3755182/files/agora103_genus.qza?download=1`
- a MICOM-formatted diet/medium: `wget -O western_diet_gut.qza https://zenodo.org/record/3755182/files/western_diet_gut.qza?download=1`
- (optional) file containing genus taxonomy translations between the tool generating the counts and the underlying taxa models used by MICOM. Column `A` shoudl have the genus as provided by your tool, and column `B` should have the AGORA genus.


## Step 0: Checking taxonomic coverage

First, run the workflow with the `--dry-run` flag to prevent snakemake from submitting any jobs.  This will run the logic to match up the genera you provided with those present in the pre-formated AGORA model.  Genera not able to matched will be written to `unmatchable_genera.csv` in your output directory. In cases of naming mismatches, create a csv file containing your and AGORA's genus names as columns `A` and `B`, and supply in the next step.

```
snakemake --config stage=tradeoff  abundances=$PWD/data.csv agora=$PWD/resources/agora103_genus.qza medium=$PWD/resources/western_diet_gut_agora.qza  --directory $PWD/results/  --dry-run
```

## Step 1: Identifying an optimal tradeoff

We will utilize Snakemake config arguments to first submit jobs calculating the optimal tradeoff parameter for each sample.  The config can be modified to set the appropriate taxa abundance cutoff.

```
snakemake --config stage=tradeoff cutoff=0.05  abundances=$PWD/data.csv agora=$PWD/resources/agora103_genus.qza medium=$PWD/resources/western_diet_gut_agora.qza  --directory $PWD/results/
```


## Step 2: Grow
Having decided on the tradeoff value based on the visualizations in `results/tradeoff.html`, you can now submit the growth jobs! This is done by changing the `stage` config argument, setting the `tradeoff` argument, and re-using the same output `--directory`.

```
snakemake --config stage=grow tradeoff=0.6 cutoff=0.05 abundances=$PWD/data.csv agora=$PWD/resources/agora103_genus.qza medium=$PWD/resources/western_diet_gut_agora.qza  --directory $PWD/results/
```

## Step 3: Analysis

From here, you can follow along the tutorial! Run the following in place of the `growth = grow(com, "models", medium, tradeoff=0.8, threads=2)` line:

```python
from micom.workflows import load_results
growth = load_results("results/growth.zip")
```
