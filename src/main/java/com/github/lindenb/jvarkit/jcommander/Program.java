/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.jcommander;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.Target;

@Retention(java.lang.annotation.RetentionPolicy.RUNTIME)
@Target({ ElementType.TYPE })
public @interface Program {
	public String name() default "";
	public String description() default "";
	public String deprecatedMsg() default "";
	public String[] keywords() default {};
	public String[] authors() default "Pierre Lindenbaum PhD @yokofakun";
	public int[] biostars() default {};
	/** shall we generate the Markdown Documentation ? */
	public boolean generate_doc() default true;
	/** bibliographic references for this paper, this paper was published in... */
	public String[] references() default {};
	/** creation date , if any */
	public String creationDate() default "";
	/** modification date , if any */
	public String modificationDate() default "";
	/** part of jvarkit amalgamion */
	public boolean jvarkit_amalgamion() default false;
	/**hidden from jvarkit amalgamion */
	public boolean jvarkit_hidden() default false;
	/** menu for jvarkit documentation */
	public String menu() default "Unclassfied";
	/**URL for nfcore */
	public String nfcore() default "";
}
