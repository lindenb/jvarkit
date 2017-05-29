# Biostar3654

show blast alignment with annotations


## Usage

```
Usage: biostar3654 [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -L, --length
      Fasta Line kength
      Default: 50
    -o, --out
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * blast
 * xml
 * annotation



## See also in Biostars

 * [https://www.biostars.org/p/3654](https://www.biostars.org/p/3654)


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
$ make biostar3654
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar3654.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar3654.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar3654** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ cat ~/jeter.blastn.xml 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|14971104|gb|AF338247.1|</Hit_id>
  <Hit_def>Human rotavirus A strain M clone M1 NSP3 genes, complete cds</Hit_def>
  <Hit_accession>AF338247</Hit_accession>
  <Hit_len>2032</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
```

```
$ java -jar dist/biostar3654.jar ~/jeter.blastn.xml 2> /dev/null  | cut -c-${COLUMNS} 

QUERY: No definition line
       ID:Query_186611 Len:980
>Human rotavirus A strain M clone M1 NSP3 genes, complete cds
 AF338247
 id:gi|14971104|gb|AF338247.1| len:2032

   e-value:0 gap:0 bitScore:1764.98

QUERY 000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ################################################## source organi
                ##################################                 5'UTR
                                                  ################ CDS codon_sta

QUERY 000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ################################################## source organi
                ################################################## CDS codon_sta
(...)

```


## Example

```
$ cat ~/jeter.blastn.xml 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|14971104|gb|AF338247.1|</Hit_id>
  <Hit_def>Human rotavirus A strain M clone M1 NSP3 genes, complete cds</Hit_def>
  <Hit_accession>AF338247</Hit_accession>
  <Hit_len>2032</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
```

```
$ java -jar dist/biostar3654.jar ~/jeter.blastn.xml 2> /dev/null  | cut -c-${COLUMNS} 

QUERY: No definition line
       ID:Query_186611 Len:980
>Human rotavirus A strain M clone M1 NSP3 genes, complete cds
 AF338247
 id:gi|14971104|gb|AF338247.1| len:2032

   e-value:0 gap:0 bitScore:1764.98

QUERY 000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ################################################## source organi
                ##################################                 5'UTR
                                                  ################ CDS codon_sta

QUERY 000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ################################################## source organi
                ################################################## CDS codon_sta
(...)
```


