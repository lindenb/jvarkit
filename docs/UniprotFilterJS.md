# UniprotFilterJS

Filters Uniprot DUMP+ XML with a javascript  (java rhino) expression. Context contain 'entry' an uniprot entry and 'index', the index in the XML file.


## Usage

```
Usage: uniprotfilterjs [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --version
      print version and exit
    -e
       (js expression). Optional.
    -f
       (js file). Optional.

```


## Keywords

 * unitprot
 * javascript
 * xjc
 * xml


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
$ make uniprotfilterjs
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/UniprotFilterJS.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/UniprotFilterJS.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **uniprotfilterjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example
the following script get the human (id=9606) uniprot entries having an id in ensembl:

```javascript
function accept(e)
	{
	var ok=0,i;
	// check organism is human 
	if(e.getOrganism()==null) return false;
	var L= e.getOrganism().getDbReference();
	if(L==null) return false;
	for(i=0;i<L.size();++i)
		{
		if(L.get(i).getId()=="9606") {ok=1;break;}
		}
	if(ok==0) return false;
	ok=0;
	L= e.getDbReference();
	if(L==null) return false;
	for(i=0;i<L.size();++i)
		{
		if(L.get(i).getType()=="Ensembl") {ok=1;break;}
		}
	return ok==1;
	}
accept(entry);
```


```bash
$   curl -skL "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz" | gunzip -c |\
java -jar dist/uniprotfilterjs.jar  -f filter.js > output.xml
```




