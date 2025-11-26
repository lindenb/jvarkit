# Json2Xml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

convert json to xml


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar json2xml  [options] Files

Usage: json2xml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --namespace, --ns, -n
      xml namespace (empty= no ns)
      Default: http://www.ibm.com/xmlns/prod/2009/jsonx
    --omit-xml-declaration
      Omit XML declaration
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * json
 * xml



## Creation Date

20251226

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/json2xml/Json2Xml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/json2xml/Json2Xml.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/json2xml/Json2XmlTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/json2xml/Json2XmlTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **json2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ echo '[null,true,{"hello":1}]' | java -jar dist/jvarkit.jar json2xml | xmllint --format -
<?xml version="1.0" encoding="UTF-8"?>
<array xmlns="http://www.ibm.com/xmlns/prod/2009/jsonx">
  <null/>
  <boolean>true</boolean>
  <object>
    <number name="hello">1</number>
  </object>
</array>

$ echo '[null,true,{"hello":1}]' | java -jar dist/jvarkit.jar json2xml --ns "" --omit-xml-declaration 
<array><null/><boolean>true</boolean><object><number name="hello">1</number></object></array>

```




