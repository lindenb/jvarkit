Notes about packaging.

See https://docs.conda.io/projects/conda-build/en/latest/user-guide/tutorials/build-pkgs.html

(supprimer `http_proxy` && `https_proxy`)

```
conda create -n ANACONDA conda-build
```

```
conda activate ANACONDA
conda-build jvarkit
conda install --use-local  jvarkit
```

```
anaconda login
anaconda upload /home/lindenb/anaconda2/envs/ANACONDA/conda-bld/linux-64/jvarkit-v20200206-0.tar.bz2
```

```
conda remove --use-local jvarkit
```
