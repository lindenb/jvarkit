# SwingPLinkSelectCluster

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Swing-based Plink/MDS sample selector


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar swingplinkselectcluster  [options] Files

Usage: swingplinkselectcluster [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --samples-to-groups, -m
      a tab delimited file with two columns: (sample-name)(TAB)(group-name)
    --version
      print version and exit

```


## Keywords

 * plink
 * sample
 * swing



## Creation Date

20231123

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/plink/SwingPLinkSelectCluster.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/plink/SwingPLinkSelectCluster.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **swingplinkselectcluster** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


GUI selecting the samples of a MDS file generated with plink.
At the end a file with fid/iid/keep-status is saved.



## Input

input can be generated with plink:

```

## Example:

```nextflow
	plink --bcf '${genome_bcf}' \\
		--double-id \\
		--read-genome '${genome_plink}' \\
		--mds-plot ${num_components} \\
		--cluster \\
		--out TMP/cluster

```

input file must have this header:

```
FID	IID	SOL	C1	C2	C3
```

## Usage

```
java -jar dist/jvarkit.jar swingplinkselectcluster file.mds  --samples-to-groups sample2group.tsv
```


