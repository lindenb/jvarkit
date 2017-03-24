/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum PhD

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2016 creation

*/
package com.github.lindenb.jvarkit.util.command;

import java.io.File;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class Arguments {
private static boolean print_hidden_options=false;

private final List<ArgumentDef> argumentdefs = new ArrayList<>();
private Predicate<ArgumentDef> argumentFilter=null;
private static final Arguments INSTANCE = new Arguments();

public static Arguments getInstance() {
	return INSTANCE;
	}

private Arguments() {
	register(Arguments.class);
	}

public ArgumentDef findArgumentByOpt(final String opt) {
	if(opt==null) return null;
	for(final ArgumentDef ad: this.argumentdefs) {
		if(opt.equals(ad.getOpt())) return ad;
	}
	return null;
}

public ArgumentDef findArgumentDefByLongOpt(final String longopt) {
	if(longopt==null) return null;
	for(final ArgumentDef ad: this.argumentdefs) {
		if(ad.hasLongOpt() && longopt.equals(ad.getLongOpt())) return ad;
	}
	return null;
}

public Arguments register(final Object obj) {
	if(obj!=null) {
		register(obj,obj.getClass());
		}
	return this;
	}

public Arguments register(final Class<?> clazz) {
	register(null,clazz);
	return this;
	}

public Arguments filter(Predicate<ArgumentDef> filter) {
	this.argumentFilter  = filter;
	int i=0;
	while(this.argumentFilter!=null && i< this.argumentdefs.size()) {
		if( !this.argumentFilter.test(this.argumentdefs.get(i))) {
			this.argumentdefs.remove(i);
			}
		else
			{
			++i;
			}
		}
	return this;
	}


private void register(final Object obj,final Class<?> clazz) {
	/** loop over the field of the class */
	for(final Field field: clazz.getDeclaredFields()) {
		//if no object, we're only looking a static fields
		if(obj==null && !java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		//if real object, we're not looking a static fields
		if(obj!=null && java.lang.reflect.Modifier.isStatic(field.getModifiers())) continue;
		//retrieve argument
		final AbstractArgument argument  = field.getAnnotation(AbstractArgument.class);
		//no argument
		if(argument==null) continue;
		//create new argument def
		final ArgumentDef def = new  ArgumentDef(obj, field,argument);
		//check we can modify it
		if(!def.isMultiple() && java.lang.reflect.Modifier.isFinal(field.getModifiers())) {
			throw new RuntimeException("not a collection and flagged final : "+def);
			}
		//discard if a argument filter is set
		if( argumentFilter!=null && !argumentFilter.test(def)) continue;
		//discard if already registerd
		if( argumentdefs.contains(def) ) {
			System.err.println("Argument already registered: "+def);
			continue;
			}
		//add to list of arguments
		this.argumentdefs.add(def);
		}
	//observe parent class
	final Class<?> parentClass = clazz.getSuperclass();
	if(parentClass==null || parentClass == Object.class) return;
	register(obj,parentClass);
	
	}


CommandLine parseOptions(final String argcargv[]) throws ParseException {
	final Options options = new Options();
	options.addOption(Option.builder("h").longOpt("help").hasArg(false).desc("display help.").build());
	for(final ArgumentDef def: this.argumentdefs) {
		options.addOption(def.getOption());
		}
	final DefaultParser parser = new DefaultParser();
	final CommandLine cli = parser.parse(options, argcargv);
	
	for(final ArgumentDef def: this.argumentdefs) {
		if(cli.hasOption(def.getOpt())) {
			def.visited=true;
			}
		else if(def.isBoolean()) {
			}
		else 
			{
			def.parseValue(cli.getOptionValue(def.getOpt()));
			}
		}
	
	return cli;
	}

public void printHelp(final PrintWriter w) {
w.print(this.argumentdefs);
w.flush();
}


private static class Test {
@AbstractArgument
boolean f1;
@AbstractArgument
Boolean f2;
@AbstractArgument
Collection<File> f3 = new ArrayList<>();
}

public static void main(String[] args)
	{
	args=new String[]{};
	Arguments arguments = new Arguments().register(new Test());
	arguments.printHelp(new PrintWriter(System.out));
	for(ArgumentDef tt:arguments.argumentdefs)
		System.err.println(tt);
	}
}