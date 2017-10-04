# BlastMapAnnotations

Maps uniprot/genbank annotations on a blast result.


## Usage

```
Usage: blastmapannots [options] Files
  Options:
    --exclude
      Exclude uniprot/feature/type of genbank/feature/key.
      Default: []
  * -u, -g, --genbank, --uniprot
      XML sequence file Genbank.xml or uniprot.xml.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --include
      Restrict to uniprot/feature/type of genbank/feature/key.
      Default: []
    --version
      print version and exit
  * -APC
      append the sequence accession before the feature name.
      Default: false

```


## Keywords

 * blast
 * annotation
 * genbank
 * uniprot


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
$ make blastmapannots
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blastmapannots/BlastMapAnnotations.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blastmapannots/BlastMapAnnotations.java)


<details>
<summary>Git History</summary>

```
Wed May 17 14:09:36 2017 +0200 ; fix typo bioalcidae ; https://github.com/lindenb/jvarkit/commit/9db2344e7ce840df02c5a7b4e2a91d6f1a5f2e8d
Tue Apr 18 17:44:13 2017 +0200 ; javacc + samfilter ; https://github.com/lindenb/jvarkit/commit/695e7cb606ba96feeeabbd2a359aacd38cf36ae0
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Fri Jul 12 11:26:27 2013 +0200 ; vcf view in gui ; https://github.com/lindenb/jvarkit/commit/b287af69d3b9bfd3a8866231aeedc5b491d314d4
Wed Jul 10 14:11:22 2013 +0200 ; fixed error in blast-map annot / genbank ; https://github.com/lindenb/jvarkit/commit/d12459b96f62d53b35263a40334938bf8000f04d
Wed Jul 10 12:34:59 2013 +0200 ; build.dtd and  fixed error in blast-map-annot ; https://github.com/lindenb/jvarkit/commit/f1b5f928840df4c894fdf8a236e4dfabf064db2c
Tue Jul 9 22:32:27 2013 +0200 ; map blast annotations ; https://github.com/lindenb/jvarkit/commit/daf4fc237904ef2cf68f6588b93d5aded76905c5
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **blastmapannots** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

Download uniprot P04514  ( Rotavirus Non-structural protein 3 )  as <b>XML</b>
```bash
$ curl -o P04514.xml "http://www.uniprot.org/uniprot/P04514.xml"
```
Download the same P04514 as <b>fasta</b>
```bash
$ curl -o P04514.fasta "http://www.uniprot.org/uniprot/P04514.fasta"
```

<b>TblastN</b> P04514.fasta vs a RNA of NSP3 in genbank http://www.ncbi.nlm.nih.gov/nuccore/AY065842.1 and save the ouput as XML:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>tblastn</BlastOutput_program>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|18139606|gb|AY065842.1|</Hit_id>
  <Hit_def>Rhesus rotavirus nonstructural protein 3 (NSP3) gene, complete cds</Hit_def>
  <Hit_accession>AY065842</Hit_accession>
  <Hit_len>1078</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_bit-score>546.584</Hsp_bit-score>
      <Hsp_score>1407</Hsp_score>
      <Hsp_evalue>0</Hsp_evalue>
      <Hsp_query-from>1</Hsp_query-from>
      <Hsp_query-to>313</Hsp_query-to>
      <Hsp_hit-from>26</Hsp_hit-from>
      <Hsp_hit-to>964</Hsp_hit-to> <Hsp_qseq>MLKMESTQQMASSIINTSFEAAVVAATSTLELMGIQYDYNEIYTRVKSKFDYVMDDSGVKNNLLGKAATIDQALNGKFGSVMRNKNWMTDSRTVAKLDEDVNKLRMMLSSKGIDQKMRVLNACFSVKRIPGKSSSVIKCTRLMKDKIERGAVEVDDSFVEEKMEVDTVDWKSRYDQLERRFESLKQRVNEKYTTWVQKAKKVNENMYSLQNVISQQQNQIADLQNYCSKLEADLQNKVGSLVSSVEWYLKSMELPDEVKTDIEQQLNSIDTISPINAIDDLEILIRNLIHDYDRTFLMFKGLLRQCNYEYAYE</Hsp_qseq>
      <Hsp_hseq>MLKMESTQQMASSIINSSFEAAVVAATSTLELMGIQYDYNEVYTRVKSKFDLVMDDSGVKNNLIGKAITIDQALNGKFSSAIRNRNWMTDSRTVAKLDEDVNKLRIMLSSKGIDQKMRVLNACFSVKRIPGKSSSIVKCTRLMKDKLERGEVEVDDSFVEEKMEVDTIDWKSRYEQLEKRFESLKHRVNEKYNHWVLKARKVNENMNSLQNVISQQQAHINELQMYNNKLERDLQSKIGSVVSSIEWYLRSMELSDDVKSDIEQQLNSIDQLNPVNAIDDFESILRNLISDYDRLFIMFKGLLQQCNYTYTYE</Hsp_hseq>
      <Hsp_midline>MLKMESTQQMASSIIN SFEAAVVAATSTLELMGIQYDYNE YTRVKSKFD VMDDSGVKNNL GKA TIDQALNGKF S  RN NWMTDSRTVAKLDEDVNKLR MLSSKGIDQKMRVLNACFSVKRIPGKSSS  KCTRLMKDK ERG VEVDDSFVEEKMEVDT DWKSRY QLE RFESLK RVNEKY  WV KA KVNENM SLQNVISQQQ  I  LQ Y  KLE DLQ K GS VSS EWYL SMEL D VK DIEQQLNSID   P NAIDD E   RNLI DYDR F MFKGLL QCNY Y YE</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
(...)
</Iteration>
</BlastOutput_iterations>
</BlastOutput>
```
Now produce a BED file with this blast result to map the features of P04514 to AY065842.

```bash
$ java -jar dist/blastmapannots.jar I=P04514.xml B=blast.xml

AY065842	25	961	Non-structural_protein_3	943	+	25961	255,255,255	1	936	25
AY065842	34	469	RNA-binding	970	+	34	469	255,255,255	1	435	34
AY065842	472	640	Dimerization	947	+	472	640	255,255,255	1	168	472
AY065842	532	724	Interaction_with_ZC3H7B	917	+	532	724	255,255,255	1	192	532
AY065842	646	961	Interaction_with_EIF4G1	905	+	646	961	255,255,255	1	315	646
AY065842	520	733	coiled-coil_region	916	+	520	733	255,255,255	1	213	520
```

