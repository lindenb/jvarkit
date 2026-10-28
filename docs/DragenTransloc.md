# DragenTransloc

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Finds dragen translocations


## Usage

```
Usage: java -jar dist/dragentransloc.jar  [options] Files
Usage: dragentransloc [options] Files
  Options:
    --accept-filtered
      Accept FILTERed variant
      Default: false
    --bcf-output
      If this program writes a VCF to a file, The format is first guessed from 
      the file suffix. Otherwise, force BCF output. The current supported BCF 
      version is : 2.1 which is not compatible with bcftools/htslib (last 
      checked 2019-11-15)
      Default: false
    --blacklist-bed, --exclude-bed
      Do not use the variant falling in that region.
    --control, --controls
      another vcf or a text file with the '.list' suffix containing the path 
      to the controls. Any translocation will be removed if it's found in a 
      control 
    --debug
      debug
      Default: false
    --disable-discordant-contig, --disable-same-contig, --disable-same-chromosome
      accept BND where both break-end are on the same chromosome.
      Default: false
    --distance
      two BND on the same chromosome at a distance of 'x' are considered the 
      same. A distance specified as a positive integer.Commas are removed. The 
      following suffixes are interpreted : b,bp,k,kb,m,mb,g,gb
      Default: 10
    --generate-vcf-md5
      Generate MD5 checksum for VCF output.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --min-dp, --min-depth
      min value for INFO/BND_DEPTH and INFO/MATE_BND_DEPTH. -1 = ignore
      Default: -1
    -o, --out
      Output file. Optional . Default: stdout
    --skip-imprecise
      Skip variants with the INFO/IMPRECISE flag.
      Default: false
    --two-junctions-factor
      When looking for the pair (ID/MATEID) in the main vcf, there should be 
      another pair of  (ID/MATEID) for the other side of the breakend, mapping 
      on the same chromosome and at a (--distance)*'x'.Bute one of the two 
      junctions may be removed when found in a control VCF. Ignored if <1.
      Default: -1
    --version
      print version and exit

```


## Keywords

 * vcf
 * dragen
 * inversion
 * bnd
 * sv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone --recurse-submodules "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew dragentransloc
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20261028

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dragentransloc/DragenTransloc.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/dragentransloc/DragenTransloc.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/dragentransloc/DragenTranslocTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/dragentransloc/DragenTranslocTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **dragentransloc** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Motivation 

find translocations using Dragen , filter out those found in another set of VCFS (controls)


e.g. true translocation found in dragen:

```
chrA  56120849  DRAGEN:BND:33583:0:1:0:12:0:1  G  G[chrB:174700[    999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:0:12:0:0;
chrA  56120857  DRAGEN:BND:33583:0:1:1:7:0:1   G  ]chrB:174678]G    999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:1:7:0:0;
chrB  174678    DRAGEN:BND:33583:0:1:1:7:0:0   T  T[chrA:56120857[  999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:1:7:0:1;
chrB  174700    DRAGEN:BND:33583:0:1:0:12:0:0  C  ]chrA:56120849]C  999  PASS  SVTYPE=BND;MATEID=DRAGEN:BND:33583:0:1:0:12:0:1;
```

# Example

```
mkdir TMP

find /path -type f -name "*.sv.vcf.gz" > controls.list

java -jar jvarkit.jar dragentransloc \\
	--controls controls.list \\
	--distance 100 \\
	"in.vcf" > out.vcf

```

