# ApplyVelocity

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Execute apache velocity macros 


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar applyvelocity  [options] Files

Usage: applyvelocity [options] Files
  Options:
    -R, --dict-file
      Load DICT file. Syntax -J KEY file.dict.
      Default: []
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -J, --json-file
      Load JSON file. Syntax -J KEY file.json.
      Default: []
    -m, --manifest
      TSV file with the following required header: key/value/type
    -o, --output
      Output file. Optional . Default: stdout
    -V, --vcf-file
      Load VCF file. Syntax -J KEY file.vcf.
      Default: []
    --version
      print version and exit

```


## Keywords

 * velocity
 * json



## Creation Date

20241023

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/velocity/ApplyVelocity.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/velocity/ApplyVelocity.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **applyvelocity** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



text utility using apache velocity templates to generate text.

## Syntax

```
cat template.vm | java -jar jvarkit.jar applyvelocity -m manifest.tsv > out.txt
cat template.vm | java -jar jvarkit.jar applyvelocity -J mycontext file.json > out.txt
java -jar jvarkit.jar applyvelocity -m manifest.tsv template1.vm template2.vm > out.txt
```

## Manifest

TSV file with 3 columns

- key :the name of the context injected in velocity
- type : type of value: 'int', 'long', 'float', 'double', 'string', 'boolean' , 'json' (might be a json string (value starts with '{' or '[' ) or a path to a json file ) , 'dict' htsjdk dictionary,
  'vcf' the variants in a vcf file. 'class' a java class, 'instance-of' instance of given java class
- value: the value 


##Example





