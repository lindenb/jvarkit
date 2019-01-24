# WesCnvTView

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

SVG visualization of bam DEPTH for multiple regions in a terminal


## Usage

```
Usage: wescnvtview [options] Files
  Options:
    -l, -B, --bams
      The Bam file(s) to be displayed. If there is only one file which ends 
      with '.list' it is interpreted as a file containing a list of paths
      Default: []
    -cap, --cap
      Cap coverage to this value. Negative=don't set any limit
      Default: -1
    -x, --extend
      Extend intervals by factor 'x'
      Default: 1.0
    --filter
      A filter expression. Reads matching the expression will be filtered-out. 
      Empty String means 'filter out nothing/Accept all'. See https://github.com/lindenb/jvarkit/blob/master/src/main/resources/javacc/com/github/lindenb/jvarkit/util/bio/samfilter/SamFilterParser.jj 
      for a complete syntax.
      Default: Accept All/ Filter out nothing
    -F, --format
      input format. INTERVALS is a string 'contig:start-end'.
      Default: INTERVALS
      Possible Values: [VCF, BED, INTERVALS]
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
    -o, --output
      Output file. Optional . Default: stdout
    -p, -percentile, --percentile
      How to compute the percentil of a region
      Default: MEDIAN
      Possible Values: [MIN, MAX, MEDIAN, AVERAGE, RANDOM, SUM]
    -P, --plain
      Plain output (not color)
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
 * svg
 * cnv


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew wescnvtview
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/WesCnvTView.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/WesCnvTView.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/WesCnvTViewTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/WesCnvTViewTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **wescnvtview** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

Input is a set of regions to observe. It can be a 

   * BED file
   * VCF File with SV
   * Some intervals 'contig:start-end'

Input can be read from stdin

## Example

```
find src/test/resources/ -type f -name "S*.bam" > bam.list

$  java -jar dist/wescnvtview.jar -l bam.list  -P "RF01:100-200" 

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
java -jar dist/wescnvtview.jar --bams bam.list -P -F BED jeter.txt   |\
	csplit -b  '%05d.txt' -f cnv. -n 5 -s -z - '/^>>>/' '{*}' 
```

`-s`: "do not print counts of output file sizes". 

`-z`: "remove empty output files". 

`-n`: "use specified number of digits instead of 2". 

`-b`: "use sprintf FORMAT instead of %02d". 

`-f`: "prefix". 

## Note to self: view in less/more

```
java -jar dist/wescnvtview.jar --bams bam.list -P -F BED jeter.txt   | less -r
```

## Screenshot

https://twitter.com/yokofakun/status/1053185975923470337

![https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg](https://pbs.twimg.com/media/Dp2rDfsWoAEQEAI.jpg)

https://twitter.com/yokofakun/status/1053204927202369536

![https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg](https://pbs.twimg.com/media/Dp28R1VWwAA7frV.jpg)

https://twitter.com/yokofakun/status/1057627022665502721

![https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg](https://pbs.twimg.com/media/Dq1x60NVAAUhGPG.jpg)

