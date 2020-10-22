# VcfBraiding

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

visualization for variants and attributes using https://visdunneright.github.io/sequence_braiding/docs/ .


## Usage

```
Usage: vcfbraiding [options] Files
  Options:
    --alleles
      show alleles in header.
      Default: false
    -B, --base
      Base URL for code 
      :https://visdunneright.github.io/sequence_braiding/docs/ 
      Default: https://visdunneright.github.io/sequence_braiding/docs/
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --hom-ref, -hr
      remove sample that are all HOM_REF for all variants.
      Default: false
    --id
      id for svg element.
      Default: vcfid
    --no-call, -nc
      remove sample that are all NO_CALL for all variants.
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -T, --title
      title
      Default: <empty string>
    --version
      print version and exit
    -D
      Dynamic parameters for API 'options'. -Dkey=value . Keys are currently: 
      show_seq_names forceLevelName animate  colorbysequence width height 
      fontSize padding.
      Syntax: -Dkey=value
      Default: {}

```


## Keywords

 * vcf
 * visualization


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfbraiding
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20201021

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/braiding/VcfBraiding.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/braiding/VcfBraiding.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfbraiding** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

## Example

```
bcftools view src/test/resources/rotavirus_rf.vcf.gz "RF02" "RF03" |\
	java -jar /home/lindenb/src/jvarkit-git/dist/vcfbraiding.jar --title "Rotavirus Variants" > variants.html
```

## Screenshots

  https://twitter.com/yokofakun/status/1319221221611941889
  
  ![https://twitter.com/yokofakun/status/1319221221611941889](https://pbs.twimg.com/media/Ek7QRz9WMAApITu?format=jpg&name=large)

  https://twitter.com/yokofakun/status/1319228442043387905
  
  ![https://twitter.com/yokofakun/status/1319228442043387905](https://pbs.twimg.com/media/Ek7W9jaXEAEQ_TN?format=png&name=small)



## See also:
 
  * https://visdunneright.github.io/sequence_braiding/docs/

