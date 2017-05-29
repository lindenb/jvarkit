# VcfGeneOntology

Find the GO terms for VCF annotated with SNPEFF or VEP


## Usage

```
Usage: vcfgo [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit
  * -A
      (goa input url)
      Default: http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD
    -C
      (Go:Term) Add children to the list of go term to be filtered. Can be 
      used multiple times.
      Default: []
    -F
       if -C is used, don't remove the variant but set the filter
  * -G
      (go  url)
      Default: http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz
    -T
      INFO tag.
      Default: GOA
    -r
      remove variant if no GO term is found associated to variant
      Default: false
    -v
      inverse filter if -C is used
      Default: false

```


## Keywords

 * vcf
 * go


## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make vcfgo
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgo/VcfGeneOntology.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfgo/VcfGeneOntology.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfgo** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example


```make
INPUT=${HOME}/test.vcf.gz
VCFGO=java -jar dist/vcfgo.jar -A  gene_association.goa_human.gz -G go_daily-termdb.rdf-xml.gz 
DEPS=dist/vcfgo.jar gene_association.goa_human.gz go_daily-termdb.rdf-xml.gz 

all: $(addsuffix .vcf,$(addprefix test,1 2 3 4 5 6 7))
	$(foreach F,$^,echo -e "Output of $F:"  && cat $F; )

test1.vcf: ${DEPS}
	${VCFGO} ${INPUT} |  grep -v "#" | grep GO | cut -f 7,8 | head -n 1 > $@

test2.vcf: ${DEPS}
	${VCFGO}  -T MYGOATAG ${INPUT} |  grep -v "#" | grep MYGOATAG | cut -f 7,8 | head -n 1 > $@
	
test3.vcf: ${DEPS}
	${VCFGO}  -C GO:0007283 ${INPUT} |  grep -v "#" | grep GO  | cut -f 7,8 | head -n 1 > $@
	
test4.vcf: ${DEPS}
	${VCFGO} -C GO:0007283 -v  ${INPUT}|  grep -v "#" |  grep GO | cut -f 7,8 | head -n 1 > $@

test5.vcf: ${DEPS}
	${VCFGO}  -C GO:0007283 -F GOFILTER ${INPUT}| grep -v "#" |   grep GOFILTER | cut -f 7,8 | head -n1 > $@
	
test6.vcf: ${DEPS}
	${VCFGO} -C GO:0007283 -F GOFILTER -v  ${INPUT}| grep -v "#" |   grep GOFILTER | cut -f 7,8 | head -n1 > $@	
	
test7.vcf: ${DEPS}
	${VCFGO} -r ${INPUT} | grep -v "#" | cut -f 7,8 | head -n1 > $@
		
	
gene_association.goa_human.gz :
	curl -o $@ "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD"

go_daily-termdb.rdf-xml.gz :
	curl -o $@ "http://archive.geneontology.org/latest-termdb/$@"

dist/vcfgo.jar:
	$(MAKE) vcfgo

```
output of GNU make:

Output of test1.vcf:

```
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231;MQ=39	GT:PL:DP:GQ	1/1:35,3,0:1:4
```

Output of test2.vcf:

```
.	(...)MYGOATAG=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test3.vcf:

```
.	(...)GOA=GJA9|GO:0007154&GO:0005922&GO:0016021,MYCBP|GO:0005739&GO:0005515&GO:0006355&GO:0005813&GO:0005737&GO:0003713&GO:0005634&GO:0007283&GO:0006351;MQ=46
```

Output of test4.vcf:
````
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test5.vcf:

```
GOFILTER	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test6.vcf:

```
GOFILTER	(...)GOA=GJA9|GO:0007154&GO:0005922&GO:0016021,MYCBP|GO:0005739&GO:0005515&GO:0006355&GO:0005813&GO:0005737&GO:0003713&GO:0005634&GO:0007283&GO:0006351
```

Output of test7.vcf:

```
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

