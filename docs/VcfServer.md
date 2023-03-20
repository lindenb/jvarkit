# VcfServer

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Web Server displaying VCF file. A web interface for vcf2table


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar vcfserver  [options] Files

Usage: vcfserver [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -jexl, --jexl
      Use/Show JEXL filter instead of Javascript filter (which is not 
      filesystem-safe). 
      Default: false
    -p, --ped, --pedigree
      Optional Pedigree file:A pedigree file. tab delimited. Columns: 
      family,id,father,mother, 
      sex:(0:unknown;1|male|M:male;2|female|F:female), phenotype 
      (-9|?|.:unknown;1|affected|case:affected;0|unaffected|control:unaffected) 
    -P, --port, -port
      Server listening port
      Default: 8080
    -timeout, --timeout
      query timeout in seconds
      Default: 60
    --url
      A custom URL for a web browser. The following words will be replaced by 
      their values: ${CHROM}, ${START}, ${END}. For example for IGV that would 
      be: 'http://localhost:60151/goto?locus=${CHROM}%3A${START}-${END}' (see 
      http://software.broadinstitute.org/software/igv/book/export/html/189) 
    --version
      print version and exit

```


## Keywords

 * vcf
 * table
 * visualization
 * server
 * web



## Creation Date

20171027

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfserver/VcfServer.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfserver/VcfServer.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfserver** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Input

Input is a set of indexed VCF file (tabix or tribble) or a file containing the path to the VCFs.


## Screenshot

https://twitter.com/yokofakun/status/923870331659485184

![twitter](https://pbs.twimg.com/media/DNI-7GZX0AA41qF.jpg "Screenshot")

https://twitter.com/yokofakun/status/924307836968079361

![twitter](https://pbs.twimg.com/media/DNPNBdQWsAAF-5w.jpg "Screenshot")



## Example 

```
$ java -jar dist/jvarkit.jar vcfserver input.vcf.gz

2017-10-27 23:53:04.140:INFO::main: Logging initialized @510ms
[INFO][VcfServer]Starting com.github.lindenb.jvarkit.tools.vcfserver.VcfServer on http://localhost:8080
2017-10-27 23:53:04.223:INFO:oejs.Server:main: jetty-9.3.7.v20160115
2017-10-27 23:53:04.336:INFO:oejs.ServerConnector:main: Started ServerConnector@9a8472{HTTP/1.1,[http/1.1]}{0.0.0.0:8080}
2017-10-27 23:53:04.337:INFO:oejs.Server:main: Started @717ms

```



