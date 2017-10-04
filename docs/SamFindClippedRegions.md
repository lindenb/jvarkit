# SamFindClippedRegions




## Usage

```
Usage: samfindclippedregions [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    --version
      print version and exit
    -B
      bed file
    -c
      min size of clipped read
      Default: 20

```

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
$ make samfindclippedregions
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

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegions.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/structvar/SamFindClippedRegions.java)


<details>
<summary>Git History</summary>

```
Fri Jun 2 16:31:30 2017 +0200 ; circos / lumpy ; https://github.com/lindenb/jvarkit/commit/7bddffca3899196e568fb5e1a479300c0038f74f
Mon May 15 12:10:21 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/b4895dd40d1c34f345cd2807f7a81395ba27e8ee
Wed May 3 17:57:20 2017 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/db456cbf0b6586ea60a4fe8ea05a5af7457d5d6e
Fri Jun 17 13:56:39 2016 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/865252a44fc018f46b4280788cec65a1383dcc18
Fri Jun 5 12:42:21 2015 +0200 ; cont ; https://github.com/lindenb/jvarkit/commit/cc909f9f4ceea181bb65e4203e3fdbde176c6f2f
Tue Feb 24 16:43:03 2015 +0100 ; vcfin : code rewrittern. picky with ALT alleles. #tweet ; https://github.com/lindenb/jvarkit/commit/65ef7741539e89c7a1a1f9cca28c13d531902c96
Fri May 23 15:00:53 2014 +0200 ; cont moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/81f98e337322928b07dfcb7a4045ba2464b7afa7
Mon May 12 14:06:30 2014 +0200 ; continue moving to htsjdk ; https://github.com/lindenb/jvarkit/commit/011f098b6402da9e204026ee33f3f89d5e0e0355
Mon May 12 10:28:28 2014 +0200 ; first sed on files ; https://github.com/lindenb/jvarkit/commit/79ae202e237f53b7edb94f4326fee79b2f71b8e8
Wed Mar 5 15:25:14 2014 +0100 ; continue everything ; https://github.com/lindenb/jvarkit/commit/15fc86fe1703398069d0bba3300be31c396dc818
Sun Mar 2 19:09:18 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/ca765415946f3ed0827af0773128178bc6aa2f62
Fri Feb 28 22:22:31 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/be0529670f639154a34f5b584c7e7ecd5362b4ea
Fri Feb 28 21:07:08 2014 +0100 ; cont ; https://github.com/lindenb/jvarkit/commit/38ed6b49556e7637a733bbef2409475729bc685d
Fri Feb 28 17:50:32 2014 +0100 ; scan clipped rgn ; https://github.com/lindenb/jvarkit/commit/e7f8693e7b4736384ddccd469c7af677bb1f0491
```

</details>

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **samfindclippedregions** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)




