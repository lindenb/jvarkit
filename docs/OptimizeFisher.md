# OptimizeFisher

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Optimize fisher test on VCF using genetic algo


## Usage

```
Usage: java -jar dist/optimizefisher.jar  [options] Files
Usage: optimizefisher [options] Files
  Options:
  * --cases
      path to a file containing the cases samples
    --chiasma
      Proba chiasma
      Default: 0.01
  * --controls
      path to a file containing the controls samples
    --disable
      disable factory(ies) by name. Comma separated
      Default: <empty string>
    --duration
      format: <integer>(years|week|days|hours|minutes|seconds)
      Default: 23h
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --max-window-size
      max sliding window side
      Default: 100000
    --min-window-size
      min sliding window side
      Default: 100
    --mutation
      Proba mutation
      Default: 0.01
    --n-samples-per-generation, -nsg
      Number of samples per generation.
      Default: 5
  * --output, -o
      Output directory
    --version
      print version and exit

```


## Keywords

 * vcf
 * burden
 * fisher


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew optimizefisher
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20221013

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/OptimizeFisher.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/burden/OptimizeFisher.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **optimizefisher** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

loads a (small) VCF and memory and, using a genetic-algorithm, tries, to find the best set of parameters, the best genomic window to find the lowest FisherTest for a burden test.

## Example

```
java -jar TMP/jvarkit.jar optimizefisher --cases TMP/cases.txt --controls TMP/ctrls.txt -o TMP TMP/normalized.vcf.gz
```

