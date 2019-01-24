# BlastFilterJS

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filters a BlastOutput with a javascript expression. The script injects each <Hit> as the variable 'blasthit'. The user script should return 'true' to keep the hit.


## Usage

```
Usage: blastfilterjs [options] Files
  Options:
    -e, --expression
       (js expression). Optional.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -f, --script
       (js file). Optional.
    --version
      print version and exit

```


## Keywords

 * blast
 * js
 * javascript
 * filter


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew blastfilterjs
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastFilterJS.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/blast/BlastFilterJS.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **blastfilterjs** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Examples

Filter Hit having <Hit_len> <1500


```

$ java -jar dist/blastfilterjs.jar blastn.xml  -e 'parseInt(hit.getHitLen())<1500' 2> /dev/null |\
  xmllint --format - | grep "Hit_len"
  
   <Hit_len>1492</Hit_len>
   <Hit_len>1488</Hit_len>
   <Hit_len>1477</Hit_len>
   <Hit_len>1452</Hit_len>
   <Hit_len>1430</Hit_len>
   <Hit_len>1064</Hit_len>
   <Hit_len>1283</Hit_len>
   <Hit_len>1052</Hit_len>
   <Hit_len>1272</Hit_len>
   <Hit_len>693</Hit_len>
     

```


keep hsp having 100 > Hsp_align-len <= 200 


```

$ cat filter.js

// keep hsp having 100>= Hsp_align-len <= 200 
function rmhsps()
	{
	var hsps = hit.getHitHsps().getHsp();
	var i=0;
	while(i< hsps.size())
		{
		var hsp = hsps.get(i);
		var hsplen = parseInt(hsp.getHspAlignLen());
		
		if( hsplen < 100 || hsplen > 300 )
			{
			hsps.remove(i);
			}
		else
			{
			i++;
			}
		}
	return true;
	}
rmhsps();

```






```

$ java -jar dist/blastfilterjs.jar -f filter.js blastn.xml 2> /dev/null |\
	xmllint --format - | grep -F 'Hsp_align-len'

	 <Hsp_align-len>289</Hsp_align-len>
	 <Hsp_align-len>291</Hsp_align-len>
	 <Hsp_align-len>197</Hsp_align-len>
	 <Hsp_align-len>227</Hsp_align-len>

```





### See also


 *  






