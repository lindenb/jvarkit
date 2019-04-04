# SkipXmlElements

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filter XML elements with a javascript  (java rhino) expression. Context contain 'element' the current element. It implementsthe interface Tag described in  SkipXmlElements.class


## Usage

```
Usage: skipxmlelements [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
    -e
       (js expression). Optional.
    -f
       (js file). Optional.

```


## Keywords

 * xml
 * javascript


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew skipxmlelements
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SkipXmlElements.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SkipXmlElements.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **skipxmlelements** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```bash
$  java -jar dist-1.128/skipxmlelements.jar -e 'element.depth < 5;'  out.blastn.xml 
```
output:
```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.30+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997)
, "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>GPIPE/9606/current/all_top_level GPIPE/9606/current/rna</BlastOutput_db>
  <BlastOutput_query-ID>1113</BlastOutput_query-ID>
  <BlastOutput_query-def>gi|11051787</BlastOutput_query-def>
  <BlastOutput_query-len>728</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>2</Parameters_sc-match>
      <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
      <Parameters_gap-open>5</Parameters_gap-open>
      <Parameters_gap-extend>2</Parameters_gap-extend>
      <Parameters_filter>L;R -d repeatmasker/repeat_9606;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>1113</Iteration_query-ID>
      <Iteration_query-def>gi|11051787</Iteration_query-def>
      <Iteration_query-len>728</Iteration_query-len>
      <Iteration_hits>
</Iteration_hits>
      <Iteration_stat>
  </Iteration_stat>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>1114</Iteration_query-ID>
      <Iteration_query-def>gi|10812511</Iteration_query-def>
      <Iteration_query-len>331</Iteration_query-len>
      <Iteration_hits>
</Iteration_hits>
      <Iteration_stat>
  </Iteration_stat>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>3</Iteration_iter-num>
      <Iteration_query-ID>1115</Iteration_query-ID>
      <Iteration_query-def>gi|10811235</Iteration_query-def>
      <Iteration_query-len>629</Iteration_query-len>
      <Iteration_hits>
</Iteration_hits>
      <Iteration_stat>
  </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```
