# FastqGrep

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Grep reads names in fastq


## DEPRECATED

use picard

## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar fastqgrep  [options] Files

Usage: fastqgrep [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -md5, --md5
      write md5 file
      Default: false
    -o, --out, -R1
      Output file for R1 fastq record or interleaved output.Output file. 
      Optional . Default: stdout
    --paired
      assume input is paired end: we expect two files, or the input is assumed 
      interleaved fastq.
      Default: false
    --version
      print version and exit
    -R
      add the read
      Default: []
    -R2
      Output file for R2 fastq record. If input is paired but R2 is omitted, 
      output will be interleaved.
    -V
      invert)
      Default: false
    -f
       file containing a list of read names
    -n
      when found, remove the read from the list of names when found more that 
      'n' time (increase speed)
      Default: -1

```


## Keywords

 * fastq


## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqGrep.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/FastqGrep.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **fastqgrep** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


TODO

## Deprecation:

use picard


