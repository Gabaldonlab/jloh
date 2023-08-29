## Install

As simple as: `git clone https://github.com/Gabaldonlab/jloh.git`
And it's ready to go! But there are a few dependencies:

| Program     | Type        | Version | Links      |
|-------------|-------------|---------|------------|
| all2vcf     | Program     | 0.7.3   | [source](https://github.com/MatteoSchiavinato/all2vcf), [cite](https://github.com/MatteoSchiavinato/all2vcf) |
| bedtools    | Program     | 2.30    | [source](https://bedtools.readthedocs.io/en/latest/), [cite](https://doi.org/10.1002/0471250953.bi1112s47) |
| Biopython   | Module      | 1.79    | [source](https://biopython.org/), [cite](https://doi.org/10.1093/bioinformatics/btp163) |
| MUMmer      | Program     | 3.1     | [source](https://anaconda.org/bioconda/mummer), [cite](https://doi.org/10.1186%2Fgb-2004-5-2-r12) |
| numpy       | Module      | 1.21.4  | [source](https://numpy.org/), [cite](https://doi.org/10.1038/s41586-020-2649-2) |
| pandas      | Module      | 1.3.5   | [source](https://pandas.pydata.org/), [cite](https://doi.org/10.5281/zenodo.3509134) |
| pybedtools  | Module      | 0.8.2   | [source](https://daler.github.io/pybedtools/main.html), [cite](https://doi.org/10.1093/bioinformatics/btr539) |
| pandarallel | Module      | 1.6.1   | [source](https://pypi.org/project/pandarallel/1.6.1/), [cite](https://github.com/nalepae/pandarallel) | 
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.6.1   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |
| samtools    |  Program    | 1.13    | [source](http://www.htslib.org/), [cite](https://doi.org/10.1093/gigascience/giab008) |

Note that **pybedtools** will look for [bedtools](https://bedtools.readthedocs.io/en/latest/) in the `$PATH`, while **pysam** will look for [samtools](http://www.htslib.org/).

The installation of **MUMmer** can be easily done via conda (`conda install -c bioconda mummer`). This will place in the `$PATH` all the toolkit from the MUMmer arsenal, in particular the tools needed by JLOH to run: nucmer, delta-filter, show-snps.

The installation of **all2vcf** is very straightforward. First clone the repository:

```
git clone https://github.com/MatteoSchiavinato/all2vcf
```

Then make a symbolic link of the main `all2vcf` executable (not of the whole folder!) into your `/bin`:

```
ln -s /path/to/all2vcf /path/to/bin
```

## Run

JLOH has different tools which are used to extract and perform different operations on LOH blocks. To see all available tools simply run:

```
./jloh -h
```

### Nextflow workflow

Together with the **JLOH** tool, we provide also a [Nextflow](http://nextflow.io/) workflow that you can use to run your samples directly from raw reads to LOH blocks. All you have to do is to edit the configuration file of the workflow (\*config) and the running script (\*sh). Edit them according to your own computer / server, and then run:

`bash reads_to_LOH_blocks.sh`

This will execute the command contained in the `*.sh` script, which in turn is a `nextflow run` command that runs the `*.nf` script with the `*.config` configuration file. In case you're working on a cluster with a slurm queuing system, you can edit the `#SBATCH` lines at the beginning and then run:

`sbatch reads_to_LOH_blocks.sh`
