JVARKIT
=======

Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.
Version     : 8ebd2be2
Compilation : 20240424122145
Github      : https://github.com/lindenb/jvarkit
Issues      : https://github.com/lindenb/jvarkit/issues

## Usage

```
  java -jar jvarkit.jar [options]
```
or
```
  java -jar jvarkit.jar <command name> (other arguments)
```

## Options

 + --help show this screen
 + --help-all show all commands, including the private ones.
 + --version print version

## Compilation

### Requirements / Dependencies

* java [compiler SDK 17](https://jdk.java.net/17/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew jvarkit
```

The java jar file will be installed in the `dist` directory.

## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **jvarkit** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

