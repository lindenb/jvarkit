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
    --ncbi-api-key
      NCBI API Key see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ 
      . If undefined. Will try to read in that order: 1) A java XML property 
      file ${HOME}/.ncbi.properties and key api_key 2) the jvm property 
      "ncbi.api.key" 3) environment variable NCBI_API_KEY
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

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) and avoid OpenJdk, use the java from Oracle. Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar3654.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar3654.java)


<details>
<summary>Git History</summary>

```
Thu Nov 2 19:54:56 2017 +0100 ; added NCBI API key ; https://github.com/lindenb/jvarkit/commit/fa13648014a42cd307b25f8661385e9f62d42bea
Mon May 29 12:33:45 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/870be8e90d7e98d947f73e67ef9965f12f351846
Wed May 24 17:27:28 2017 +0200 ; lowres bam2raster & fix doc ; https://github.com/lindenb/jvarkit/commit/6edcfd661827927b541e7267195c762e916482a0
Thu May 11 16:20:27 2017 +0200 ; move to jcommander ; https://github.com/lindenb/jvarkit/commit/15b6fabdbdd7ce0d1e20ca51e1c1a9db8574a59e
Wed Apr 12 07:17:41 2017 +0200 ; biostar3654 fixed ; https://github.com/lindenb/jvarkit/commit/fe9b6b5c09e2456a16f99944d733531793b37ffa
Thu Jul 28 09:48:29 2016 +0200 ; NCBI moved API to https ; https://github.com/lindenb/jvarkit/commit/d207e023a06d2ae7afd2e05d2f1369b8a713974b
Wed Apr 22 12:21:29 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/dc65752a1d0c364957940847f8901d32106f21c7
Wed Apr 22 00:38:57 2015 +0200 ; Biostar3654.java copied from http://plindenbaum.blogspot.fr/2010/11/blastxmlannotations.html ; https://github.com/lindenb/jvarkit/commit/f57fc24c79f0a542edc85b6fe9c0900b64b0c3bb
```

</details>

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


