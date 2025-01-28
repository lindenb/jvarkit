# WGSCoveragePlotter

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Whole genome coverage plotter


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar wgscoverageplotter  [options] Files

Usage: wgscoverageplotter [options] Files
  Options:
    --clip, --cap
      Don't show coverage to be greater than 'max-depth' in the SVG file.
      Default: false
    --dimension
      Image Dimension. a dimension can be specified as '[integer]x[integer]' 
      or the word 'screen' or it can be the path to an existing 
      png,jpg,xcf,svg file.
      Default: java.awt.Dimension[width=1000,height=500]
    --disable-paired-overlap
      Count overlapping bases with mate for paired-end
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -I, --include-contig-regex
      Only keep chromosomes matching this regular expression. Ignore if blank.
      Default: <empty string>
    --mapq
      min mapping quality
      Default: 1
    -C, --max-depth
      Max depth to display. The special value '-1' will first compute the 
      average depth and the set the max depth to 2*average
      Default: 100
    --min-contig-length
      Skip chromosome with length < 'x'. A distance specified as a positive 
      integer.Commas are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb,g,gb 
      Default: 0
    -o, --output
      Output file. Optional . Default: stdout
    --partition
      When using the option --samples, use this partition Data partitioning 
      using the SAM Read Group (see 
      https://gatkforums.broadinstitute.org/gatk/discussion/6472/ ) . It can 
      be any combination of sample, library....
      Default: sample
      Possible Values: [readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center, any]
    --percentile
      How to we bin the coverage under one pixel.
      Default: median
      Possible Values: [median, average, min, max]
    --points
      Plot the coverage using points instead of areas.
      Default: false
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --samples
      Limit to those groups. See also --partition. Multiple separated with 
      commas. 
      Default: <empty string>
    -X, --skip-contig-regex
      Skip chromosomes matching this regular expression. Ignore if blank.
      Default: <empty string>
    --version
      print version and exit
    -D
      set some css style elements. '-Dkey=value'. Undocumented.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * svg
 * bam
 * depth
 * coverage



## See also in Biostars

 * [https://www.biostars.org/p/104063](https://www.biostars.org/p/104063)
 * [https://www.biostars.org/p/475162](https://www.biostars.org/p/475162)
 * [https://www.biostars.org/p/9536274](https://www.biostars.org/p/9536274)



## NF-CORE

![nfcorelogo](https://avatars.githubusercontent.com/u/35520196?s=32&v=4) This program is available in nf-core at [https://nf-co.re/modules/jvarkit_wgscoverageplotter.html](https://nf-co.re/modules/jvarkit_wgscoverageplotter.html)


## Creation Date

20201125

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotter.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotterTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2graphics/WGSCoveragePlotterTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **wgscoverageplotter** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Files

Input is an indexed BAM or CRAM file

Output is a SVG file

## Example
```
java -jar dist/jvarkit.jar wgscoverageplotter --dimension 1500x500 -C -1 --clip -R src/test/resources/rotavirus_rf.fa src/test/resources/S1.bam --include-contig-regex "RF.*" --percentile median  > ~/jeter.svg
```

## Cited in

 * Vinicius A.C & al. Comparative analyses of Theobroma cacao and T. grandiflorum mitogenomes reveal conserved gene content embedded within complex and plastic structures". Gene Volume 849 2023 .https://doi.org/10.1016/j.gene.2022.146904.
 * Varani AM, Silva SR, Lopes S, Barbosa JBF, Oliveira D, Correa MA, Moraes AP, Miranda VFO, Prosdocimi F. 2022. The complete organellar genomes of the entheogenic plant Psychotria viridis (Rubiaceae), a main component of the ayahuasca brew. PeerJ 10:e14114 https://doi.org/10.7717/peerj.14114
 * Cho, C.H., Park, S.I., Huang, TY. et al. Genome-wide signatures of adaptation to extreme environments in red algae. Nat Commun 14, 10 (2023). https://doi.org/10.1038/s41467-022-35566-x
 * Rhys T. White, Martina Jelocnik, Natalie Klukowski, Md. Hakimul Haque, Subir Sarker. The first genomic insight into Chlamydia psittaci sequence type (ST)24 from a healthy captive psittacine host in Australia demonstrates evolutionary proximity with strains from psittacine, human, and equine hosts. Veterinary Microbiology, Volume 280, 2023, 109704, ISSN 0378-1135. https://doi.org/10.1016/j.vetmic.2023.109704.
 * Rhys T. White. A discovery down under: decoding the draft genome sequence of Pantoea stewartii from Australia's Critically Endangered western ground parrot/kyloring (Pezoporus flaviventris). https://doi.org/10.1099/mgen.0.001101 . Microbial Genomics Vol9 Issue 9
 * genomic and transcriptomic insight into the evolution and breeding of turfgrass species, poa annua and bouteloua dactyloides. Thesis. Christopher W. Benson. https://etda.libraries.psu.edu/files/final_submissions/29057 **WGSCoveragePlotter**
 * Jo, J., Park, JS., Won, H. et al. The first Chromosomal-level genome assembly of Sageretia thea using Nanopore long reads and Pore-C technology. Sci Data 11, 959 (2024). https://doi.org/10.1038/s41597-024-03798-9
 * Tim P. Bean, Hannah Farley, Jennifer Nascimento-Schulze, Tim Regan,Scottish oyster mortality event and association with Vibrio aestuarianus, Aquaculture Reports, Volume 39, 2024, 102480, ISSN 2352-5134, https://doi.org/10.1016/j.aqrep.2024.102480. **WGSCoveragePlotter**
 * del Olmo V, Redondo-Rio A, Garcia AB, Limtong S, Saus E, Gabaldon T (2025) Insights into the origin, hybridisation and adaptation of Candida metapsilosis hybrid pathogens. PLoS Pathog 21(1): e1012864. https://doi.org/10.1371/journal.ppat.1012864

## Screenshot

https://twitter.com/yokofakun/status/1331898068002861056

![twitter](https://pbs.twimg.com/media/EnvaOnNW4AAkGTz?format=jpg&name=medium "Screenshot")



