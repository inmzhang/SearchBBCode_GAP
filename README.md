# SearchBBCode_GAP

Search [IBM's Bivariate Bicycle LDPC codes](http://arxiv.org/abs/2308.07915) with GAP system.

## Usage

### Install HPC-GAP

[HPC-GAP](https://github.com/gap-system/gap/wiki/Building-HPC-GAP) should be installed for paralleling the search:

```shell
chmod +x install_hpc_gap.sh

./install_hpc_gap.sh
```

### Search

A python script `search.py` is provided for searching the codes. The script will generate the GAP scripts and run it with HPC-GAP.

```shell
$ python search.py --help

usage: search.py [-h] -l L -m M [-o OUT_DIR] [-r ERLB] [--dlb DLB] [-n NIS] [-c CHUNK_SIZE] [-t TARGET_K] [-p P] [-f]

options:
  -h, --help            show this help message and exit
  -l L
  -m M
  -o OUT_DIR, --out_dir OUT_DIR
  -r ERLB, --erlb ERLB  Encoding rate lower bound
  --dlb DLB             Distance lower bound
  -n NIS, --nis NIS     Number of information sets used for distance upper bound estimation
  -c CHUNK_SIZE, --chunk_size CHUNK_SIZE
                        Chunk size for each searching subtask
  -t TARGET_K, --target_k TARGET_K
                        Target number of logical qubits for searching
  -p P                  Number of logical processors used for hpc-gap
  -f, --format_csv      Format the code coefficients in the csv file
```

For example, to search the [[72, 12, 6]] code, we can run:

```shell
python search.py -l 6 -m 6 --erlb 1/13 --dlb 5 --nis 10 -p 8 -f
```

The search results will be saved in the `out` directory by default. The following codes have been found:

```csv
As,Bs,encodingRate,parameter
x^1+x^2+y^3,x^3+y^1+y^2,1/12,"[72, 12, 6]"
x^1+x^2+y^3,x^3+y^4+y^5,1/12,"[72, 12, 6]"
x^3+y^1+y^2,x^4+x^5+y^3,1/12,"[72, 12, 6]"
x^3+y^4+y^5,x^4+x^5+y^3,1/12,"[72, 12, 6]"
```

The speed of the search heavily depends on `--erlb`/`--nis`/`--dlb` options, which will early stop the search if the the set of code parameters violate some conditions.