# CnvTView

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Text visualization of bam DEPTH for multiple regions in a terminal


## Usage

```
Usage: cnvtview [options] Files
  Options:
    -cap, --cap
      Cap coverage to this value. Negative=don't set any limit
      Default: -1
    -x, --extend
      Extending interval. The following syntaxes are supported: 1000; 1kb; 
      1,000; 30%(shrink); 150% (extend); 0.5 (shrink); 1.5 (extend)
      Default: com.github.lindenb.jvarkit.samtools.util.IntervalExtender$ExtendByFraction fraction : 1.5
    --flush
      do not wait for all bam to be scanned, do not sort, do not normalize on 
      the depth of all bams, print the figure as soon as possible.
      Default: false
    --format
      output format
      Default: plain
      Possible Values: [plain, ansi]
    -G, --genes
      A BED file containing some regions of interest that will be displayed
    -H, --height, --rows
      Terminal width per sample
      Default: 10
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -highlight, --highlight, --top
      Per default samples are sorted alphabetically.The samples in this 
      collection will be displayed on the 'top' to have an quick insight about 
      the propositus.
      Default: []
  * -r, --regions, --interval
      A source of intervals. The following suffixes are recognized: vcf, 
      vcf.gz bed, bed.gz, gtf, gff, gff.gz, gtf.gz.Otherwise it could be an 
      empty string (no interval) or a list of plain interval separated by '[ 
      \t\n;,]' 
      Default: (empty)
    --mapq
      Min mapping quality
      Default: 1
    -o, --output
      Output file. Optional . Default: stdout
    -p, -percentile, --percentile
      How to compute the percentil of a region
      Default: MEDIAN
      Possible Values: [AVERAGE, MEDIAN]
  * -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard CreateSequenceDictionary
    --stddev
      Sort output on standard deviation
      Default: false
    --version
      print version and exit
    -w, --width, --cols, -C
      Terminal width. Under linux good idea is to use the environment variable 
      ${COLUMNS} 
      Default: 80

```


## Keywords

 * bam
 * alignment
 * graphics
 * visualization
 * cnv
 * ascii
 * text


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew cnvtview
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20181018

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/CnvTView.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/CnvTView.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/CnvTViewTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/CnvTViewTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **cnvtview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

Input is a set of indexed cram/bam files or a file with the '.list' suffix containing the path to the bams

## Example

```
find src/test/resources/ -type f -name "S*.bam" > bam.list

$  java -jar dist/cnvtview.jar  -r "RF01:100-200" bam.list 

>>> RF01:100-200	 Length:101	(1)
> S1 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                %                                                      
   6.30 |           %%%%%%%%                                                    
   5.40 |           %%%%%%%%%      #           # #                              
   4.50 |       %%%%%%%%%%%%%   % ##  ####   #########   %%%%%                  
   3.60 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%% %%%%   %%%     %
   2.70 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |   %%%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 | %%%%%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S2 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                           %%%  %%     
   8.10 |                                                           %%% %%%     
   7.20 |                                                   %%  %%%%%%%%%%%     
   6.30 |                     %                             %%%%%%%%%%%%%%%%%   
   5.40 |              %%%%%%%%   #       ###           #  %%%%%%%%%%%%%%%%%%   
   4.50 |          %%%%%%%%%%%%%%##  ##  #########    ###%%%%%%%%%%%%%%%%%%%%%%%
   3.60 |          %%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   2.70 |      %%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |     %%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S3 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                           %%%  %%     
   8.10 |                                                           %%% %%%     
   7.20 |                                                   %%  %%%%%%%%%%%     
   6.30 |                     %                             %%%%%%%%%%%%%%%%%   
   5.40 |              %%%%%%%%   #       ###           #  %%%%%%%%%%%%%%%%%%   
   4.50 |          %%%%%%%%%%%%%%##  ##  #########    ###%%%%%%%%%%%%%%%%%%%%%%%
   3.60 |          %%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   2.70 |      %%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |     %%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S4 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                                                                       
   6.30 |                                    ##                                 
   5.40 |                            ##########                                 
   4.50 |                      % ################        %  %%                  
   3.60 |                      %%################# ######%%%%%%%%%%             
   2.70 |                %%%%%%%%########################%%%%%%%%%%%%%          
   1.80 |        %%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%  %    
   0.90 |       %%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%

> S5 ===========================================================================
     Pos| 1 9 18 31 44 56 69 82 95 108 125 142 160 177 194 211 228 246 263 280  
   9.00 |                                                                       
   8.10 |                                                                       
   7.20 |                                                             %%       %
   6.30 |                                                           %%%%       %
   5.40 |                                                           %%%%%     %%
   4.50 |                       %#######                 %%%%%     %%%%%% %%%%%%
   3.60 |                   % %%%########                %%%%%%%%%%%%%%%%%%%%%%%
   2.70 |                %%%%%%%%########### #          #%%%%%%%%%%%%%%%%%%%%%%%
   1.80 |               %%%%%%%%%###############        #%%%%%%%%%%%%%%%%%%%%%%%
   0.90 |    %%%%%%%%%%%%%%%%%%%%########################%%%%%%%%%%%%%%%%%%%%%%%
<<< RF01:100-200	 Length:101	(1)


```

## Note to self: Splitting the output:

```
java -jar dist/cnvtview.jar --bams bam.list -P -F BED jeter.txt   |\
	csplit -b  '%05d.txt' -f cnv. -n 5 -s -z - '/^>>>/' '{*}' 
```

`-s`: "do not print counts of output file sizes". 

`-z`: "remove empty output files". 

`-n`: "use specified number of digits instead of 2". 

`-b`: "use sprintf FORMAT instead of %02d". 

`-f`: "prefix". 

## Note to self: view in less/more

```
java -jar dist/cnvtview.jar bam.list -r jeter.txt   | less -r
```

## Screenshot

https://twitter.com/yokofakun/status/1053185975923470337

![https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg](https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg)

https://twitter.com/yokofakun/status/1053204927202369536

![https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg](https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg)

https://twitter.com/yokofakun/status/1057627022665502721

![https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg](https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg)

