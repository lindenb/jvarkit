# Bam2Xml

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

converts a BAM to XML


## Usage

```
Usage: bam2xml [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * xml


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew bam2xml
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2xml/Bam2Xml.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/bam2xml/Bam2Xml.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2xml/Bam2XmlTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/bam2xml/Bam2XmlTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **bam2xml** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




```
 samtools view samtools-1.4/examples/toy.sam |\
 	java -jar dist/bam2xml.jar  | xmllint --format -


<?xml version="1.0" encoding="UTF-8"?>
<sam>
  <header version="1.5" sort="unsorted">
    <dict size="0"/>
    <read-groups/>
    <program-records/>
  </header>
  <record id="1" flag="163" length="19" read-paired="true" proper-pair="true" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="true" first-of-pair="false" second-of-pair="true" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="7" unclipped-align-start="7" align-end="22" unclipped-align-end="22" mapq="30" mate-ref-name="ref" mate-tid="-1" mate-align-start="37" insert-size="39">
    <name>r001</name>
    <seq>TTAGATAAAGAGGATACTG</seq>
    <cigar>
      <ce op="M" length="8" read-pos="0" ref-pos="7"/>
      <ce op="I" length="4" read-pos="8"/>
      <ce op="M" length="4" read-pos="12" ref-pos="15"/>
      <ce op="D" length="1" ref-pos="19"/>
      <ce op="M" length="3" read-pos="16" ref-pos="20"/>
    </cigar>
    <attributes>
      <attribute name="XX">[S@1e81f4dc</attribute>
    </attributes>
  </record>
  <record id="2" flag="0" length="17" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="9" unclipped-align-start="8" align-end="18" unclipped-align-end="18" mapq="30">
    <name>r002</name>
    <seq>AAAAGATAAGGGATAAA</seq>
    <cigar>
      <ce op="S" length="1" read-pos="0" ref-pos="8"/>
      <ce op="I" length="2" read-pos="1"/>
      <ce op="M" length="6" read-pos="3" ref-pos="9"/>
      <ce op="P" length="1"/>
      <ce op="I" length="1" read-pos="9"/>
      <ce op="P" length="1"/>
      <ce op="I" length="1" read-pos="10"/>
      <ce op="M" length="4" read-pos="11" ref-pos="15"/>
      <ce op="I" length="2" read-pos="15"/>
    </cigar>
    <attributes/>
  </record>
  <record id="3" flag="0" length="6" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="9" unclipped-align-start="4" align-end="14" unclipped-align-end="14" mapq="30">
    <name>r003</name>
    <seq>AGCTAA</seq>
    <cigar>
      <ce op="H" length="5" ref-pos="4"/>
      <ce op="M" length="6" read-pos="0" ref-pos="9"/>
    </cigar>
    <attributes/>
  </record>
  <record id="4" flag="0" length="12" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="16" unclipped-align-start="16" align-end="40" unclipped-align-end="40" mapq="30">
    <name>r004</name>
    <seq>ATAGCTCTCAGC</seq>
    <cigar>
      <ce op="M" length="6" read-pos="0" ref-pos="16"/>
      <ce op="N" length="14" ref-pos="22"/>
      <ce op="I" length="1" read-pos="6"/>
      <ce op="M" length="5" read-pos="7" ref-pos="36"/>
    </cigar>
    <attributes/>
  </record>
  <record id="5" flag="16" length="5" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="true" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="29" unclipped-align-start="23" align-end="33" unclipped-align-end="33" mapq="30">
    <name>r003</name>
    <seq>TAGGC</seq>
    <cigar>
      <ce op="H" length="6" ref-pos="23"/>
      <ce op="M" length="5" read-pos="0" ref-pos="29"/>
    </cigar>
    <attributes/>
  </record>
  <record id="6" flag="83" length="9" read-paired="true" proper-pair="true" read-unmapped="false" mate-unmapped="false" read-reverse-strand="true" mate-reverse-strand="false" first-of-pair="true" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref" tid="-1" align-start="37" unclipped-align-start="37" align-end="45" unclipped-align-end="45" mapq="30" mate-ref-name="ref" mate-tid="-1" mate-align-start="7" insert-size="-39">
    <name>r001</name>
    <seq>CAGCGCCAT</seq>
    <cigar>
      <ce op="M" length="9" read-pos="0" ref-pos="37"/>
    </cigar>
    <attributes/>
  </record>
  <record id="7" flag="0" length="20" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="1" unclipped-align-start="1" align-end="20" unclipped-align-end="20" mapq="30">
    <name>x1</name>
    <seq>AGGTTTTATAAAACAAATAA</seq>
    <qual>????????????????????</qual>
    <cigar>
      <ce op="M" length="20" read-pos="0" ref-pos="1"/>
    </cigar>
    <attributes/>
  </record>
  <record id="8" flag="0" length="21" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="2" unclipped-align-start="2" align-end="22" unclipped-align-end="22" mapq="30">
    <name>x2</name>
    <seq>GGTTTTATAAAACAAATAATT</seq>
    <qual>?????????????????????</qual>
    <cigar>
      <ce op="M" length="21" read-pos="0" ref-pos="2"/>
    </cigar>
    <attributes/>
  </record>
  <record id="9" flag="0" length="26" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="6" unclipped-align-start="6" align-end="27" unclipped-align-end="27" mapq="30">
    <name>x3</name>
    <seq>TTATAAAACAAATAATTAAGTCTACA</seq>
    <qual>??????????????????????????</qual>
    <cigar>
      <ce op="M" length="9" read-pos="0" ref-pos="6"/>
      <ce op="I" length="4" read-pos="9"/>
      <ce op="M" length="13" read-pos="13" ref-pos="15"/>
    </cigar>
    <attributes/>
  </record>
  <record id="10" flag="0" length="25" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="10" unclipped-align-start="10" align-end="34" unclipped-align-end="34" mapq="30">
    <name>x4</name>
    <seq>CAAATAATTAAGTCTACAGAGCAAC</seq>
    <qual>?????????????????????????</qual>
    <cigar>
      <ce op="M" length="25" read-pos="0" ref-pos="10"/>
    </cigar>
    <attributes/>
  </record>
  <record id="11" flag="0" length="24" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="12" unclipped-align-start="12" align-end="35" unclipped-align-end="35" mapq="30">
    <name>x5</name>
    <seq>AATAATTAAGTCTACAGAGCAACT</seq>
    <qual>????????????????????????</qual>
    <cigar>
      <ce op="M" length="24" read-pos="0" ref-pos="12"/>
    </cigar>
    <attributes/>
  </record>
  <record id="12" flag="0" length="23" read-paired="false" proper-pair="false" read-unmapped="false" mate-unmapped="false" read-reverse-strand="false" mate-reverse-strand="false" first-of-pair="false" second-of-pair="false" not-primary-alignment="false" read-fails-vendor-quality-check="false" duplicate-read="false" supplementary-alignment="false" ref-name="ref2" tid="-1" align-start="14" unclipped-align-start="14" align-end="36" unclipped-align-end="36" mapq="30">
    <name>x6</name>
    <seq>TAATTAAGTCTACAGAGCAACTA</seq>
    <qual>???????????????????????</qual>
    <cigar>
      <ce op="M" length="23" read-pos="0" ref-pos="14"/>
    </cigar>
    <attributes/>
  </record>
</sam>

```

