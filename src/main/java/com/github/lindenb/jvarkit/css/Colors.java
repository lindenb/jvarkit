package com.github.lindenb.jvarkit.css;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.lang.CharSplitter;

/**
 * Utilities for CSS, SVG colors
 * @author lindenb
 *
 */
public class Colors {
	@SuppressWarnings("serial")
	private static final Map<String,int[]> named= new HashMap<String,int[]>() 
			{{{
				put("aliceblue",new int[]{240, 248, 255});
				put("antiquewhite",new int[]{250, 235, 215});
				put("aquamarine",new int[]{127, 255, 212});
				put("aqua",new int[]{ 0, 255, 255});
				put("azure",new int[]{240, 255, 255});
				put("beige",new int[]{245, 245, 220});
				put("bisque",new int[]{255, 228, 196});
				put("blanchedalmond",new int[]{255, 235, 205});
				put("blueviolet",new int[]{138, 43, 226});
				put("brown",new int[]{165, 42, 42});
				put("burlywood",new int[]{222, 184, 135});
				put("blue",new int[]{0,0,255});
				put("cadetblue",new int[]{ 95, 158, 160});
				put("carmine",new int[]{215, 0, 64});
				put("cerise",new int[]{222, 49, 99});
				put("chartreuse",new int[]{127, 255, 0});
				put("chocolate",new int[]{210, 105, 30});
				put("coral",new int[]{255, 127, 80});
				put("cornflowerblue",new int[]{100, 149, 237});
				put("cornsilk",new int[]{255, 248, 220});
				put("crimson",new int[]{220, 20, 60});
				put("darkblue",new int[]{ 0, 0, 139});
				put("darkcyan",new int[]{ 0, 139, 139});
				put("darkgoldenrod",new int[]{184, 134, 11});
				put("darkgray",new int[]{169, 169, 169});
				put("darkgreen",new int[]{ 0, 100, 0});
				put("darkgrey",new int[]{169, 169, 169});
				put("darkkhaki",new int[]{189, 183, 107});
				put("darkmagenta",new int[]{139, 0, 139});
				put("darkolivegreen",new int[]{ 85, 107, 47});
				put("darkorange",new int[]{255, 140, 0});
				put("darkorchid",new int[]{153, 50, 204});
				put("darkred",new int[]{139, 0, 0});
				put("darksalmon",new int[]{233, 150, 122});
				put("darkseagreen",new int[]{143, 188, 143});
				put("darkslateblue",new int[]{ 72, 61, 139});
				put("darkslategray",new int[]{ 47, 79, 79});
				put("darkslategrey",new int[]{ 47, 79, 79});
				put("darkturquoise",new int[]{ 0, 206, 209});
				put("darkviolet",new int[]{148, 0, 211});
				put("deeppink",new int[]{255, 20, 147});
				put("deepskyblue",new int[]{ 0, 191, 255});
				put("dimgray",new int[]{105, 105, 105});
				put("dimgrey",new int[]{105, 105, 105});
				put("dodgerblue",new int[]{ 30, 144, 255});
				put("firebrick",new int[]{178, 34, 34});
				put("floralwhite",new int[]{255, 250, 240});
				put("forestgreen",new int[]{ 34, 139, 34});
				put("fuchsia",new int[]{255, 0, 255});
				put("gainsboro",new int[]{220, 220, 220});
				put("ghostwhite",new int[]{248, 248, 255});
				put("goldenrod",new int[]{218, 165, 32});
				put("gold",new int[]{255, 215, 0});
				put("gray",new int[]{128, 128, 128});
				put("green",new int[]{0,128,0});
				put("greenyellow",new int[]{173, 255, 47});
				put("grey",new int[]{128, 128, 128});
				put("honeydew",new int[]{240, 255, 240});
				put("hotpink",new int[]{255, 105, 180});
				put("indianred",new int[]{205, 92, 92});
				put("indigo",new int[]{ 75, 0, 130});
				put("ivory",new int[]{255, 255, 240});
				put("khaki",new int[]{240, 230, 140});
				put("lavenderblush",new int[]{255, 240, 245});
				put("lavender",new int[]{230, 230, 250});
				put("lawngreen",new int[]{124, 252, 0});
				put("lemonchiffon",new int[]{255, 250, 205});
				put("lightblue",new int[]{173, 216, 230});
				put("lightcoral",new int[]{240, 128, 128});
				put("lightcyan",new int[]{224, 255, 255});
				put("lightgoldenrodyellow",new int[]{250, 250, 210});
				put("lightgray",new int[]{211, 211, 211});
				put("lightgreen",new int[]{144, 238, 144});
				put("lightgrey",new int[]{211, 211, 211});
				put("lightpink",new int[]{255, 182, 193});
				put("lightsalmon",new int[]{255, 160, 122});
				put("lightseagreen",new int[]{ 32, 178, 170});
				put("lightskyblue",new int[]{135, 206, 250});
				put("lightslategray",new int[]{119, 136, 153});
				put("lightslategrey",new int[]{119, 136, 153});
				put("lightsteelblue",new int[]{176, 196, 222});
				put("lightyellow",new int[]{255, 255, 224});
				put("limegreen",new int[]{ 50, 205, 50});
				put("lime",new int[]{ 0, 255, 0});
				put("linen",new int[]{250, 240, 230});
				put("magenta",new int[]{255, 0, 255});
				put("maroon",new int[]{128, 0, 0});
				put("mediumaquamarine",new int[]{102, 205, 170});
				put("mediumblue",new int[]{ 0, 0, 205});
				put("mediumorchid",new int[]{186, 85, 211});
				put("mediumpurple",new int[]{147, 112, 219});
				put("mediumseagreen",new int[]{ 60, 179, 113});
				put("mediumslateblue",new int[]{123, 104, 238});
				put("mediumspringgreen",new int[]{ 0, 250, 154});
				put("mediumturquoise",new int[]{ 72, 209, 204});
				put("mediumvioletred",new int[]{199, 21, 133});
				put("midnightblue",new int[]{ 25, 25, 112});
				put("mintcream",new int[]{245, 255, 250});
				put("mistyrose",new int[]{255, 228, 225});
				put("moccasin",new int[]{255, 228, 181});
				put("navajowhite",new int[]{255, 222, 173});
				put("navy",new int[]{ 0, 0, 128});
				put("oldlace",new int[]{253, 245, 230});
				put("olivedrab",new int[]{107, 142, 35});
				put("olive",new int[]{128, 128, 0});
				put("orangered",new int[]{255, 69, 0});
				put("orchid",new int[]{218, 112, 214});
				put("palegoldenrod",new int[]{238, 232, 170});
				put("palegreen",new int[]{152, 251, 152});
				put("paleturquoise",new int[]{175, 238, 238});
				put("palevioletred",new int[]{219, 112, 147});
				put("papayawhip",new int[]{255, 239, 213});
				put("peachpuff",new int[]{255, 218, 185});
				put("peru",new int[]{205, 133, 63});
				put("plum",new int[]{221, 160, 221});
				put("powderblue",new int[]{176, 224, 230});
				put("purple",new int[]{128, 0, 128});
				put("red",new int[]{255,0,0});
				put("rosybrown",new int[]{188, 143, 143});
				put("royalblue",new int[]{ 65, 105, 225});
				put("saddlebrown",new int[]{139, 69, 19});
				put("salmon",new int[]{250, 128, 114});
				put("sandybrown",new int[]{244, 164, 96});
				put("seagreen",new int[]{ 46, 139, 87});
				put("seashell",new int[]{255, 245, 238});
				put("sienna",new int[]{160, 82, 45});
				put("silver",new int[]{192, 192, 192});
				put("skyblue",new int[]{135, 206, 235});
				put("slateblue",new int[]{106, 90, 205});
				put("slategray",new int[]{112, 128, 144});
				put("slategrey",new int[]{112, 128, 144});
				put("snow",new int[]{255, 250, 250});
				put("springgreen",new int[]{ 0, 255, 127});
				put("steelblue",new int[]{ 70, 130, 180});
				put("tan",new int[]{210, 180, 140});
				put("teal",new int[]{ 0, 128, 128});
				put("thistle",new int[]{216, 191, 216});
				put("tomato",new int[]{255, 99, 71});
				put("turquoise",new int[]{ 64, 224, 208});
				put("violet",new int[]{238, 130, 238});
				put("wheat",new int[]{245, 222, 179});
				put("whitesmoke",new int[]{245, 245, 245});
				put("yellow",new int[]{255,255,0});
				put("yellowgreen",new int[]{154, 205, 50});
			}}};
	
	public static Set<String> getNames() {
		return Collections.unmodifiableSet(named.keySet());
		}
	
	public static String toRGB(int[] array) {
		switch(array.length) {
			case 1: return toRGB(array[0],array[0],array[0]);
			case 3: return toRGB(array[0],array[1],array[2]);
			case 4: return toRGB(array[0],array[1],array[2],array[3]);
			default: throw new IllegalArgumentException("bad number of color components");
			}
		}

			
	public static String toRGB(int r, int g,int b) {
		return "rgb("+r+","+g+","+b+")";
		}
	
	public static String toRGB(int r, int g,int b,int a) {
		if(a==255) return toRGB(r,g,b);
		return "rgb("+r+","+g+","+b+","+a+")";
		}
	
	
	
	/** add alpha component if missing */
	private static int[] alpha(int c[]) {
		if(c.length==4) return c;
		if(c.length!=3) throw new IllegalArgumentException();
		int[] array = new int[4];
		System.arraycopy(c, 0, array, 0, 3);
		array[3]=255;
		return array;
		}

    /**
     * parse Colors as a string.
     * String can be a nominal SVG value e.g. "blue","red"... or
     * a RGB definition such as "rgb(100,200,300)" or "#ab12cc"
     * @param c the color as a string
     * @return the Color as an array of int
     */
    public static int[] componentsFor(String c)
        {
        c=c.trim().toLowerCase();
        if(c.startsWith("#"))
            {
        	int value = 0xff000000 |  Integer.valueOf(c.substring(1),16).intValue();
            return new int[] {
            		(value >> 16) & 0xFF,
            		(value >> 8) & 0xFF,
            		(value >> 0) & 0xFF,
            		(value >> 24) & 0xFF
            	};
            }
        else if((c.startsWith("rgb(") || c.startsWith("rgba(")) && c.endsWith(")"))
            {
        	final int par = c.indexOf('(');
        	final String  tokens[]= CharSplitter.COMMA.split(c.substring(par+1, c.length()-1).trim());
            if(tokens.length<3)
                    {
                    return null;
                    }
            final int r= java.lang.Integer.parseInt(tokens[0]);
            final int g= java.lang.Integer.parseInt(tokens[1]);
            final int b= java.lang.Integer.parseInt(tokens[2]);
            if(tokens.length==3) {
            	return new int[] {r,g,b};
            	}
            else
            	{
                final int a= java.lang.Integer.parseInt(tokens[3]);
                return new int[] {r,g,b,a};
            	}
            }
        int[] array = named.getOrDefault(c, null);
        if(array!=null) return array;
        throw new IllegalArgumentException("Illegal Colors:"+c);
        }
	
	public static String shadeOf(
			final float f,
			final String named1,
			final String named2
			)
		{
		final int[]  c1 = alpha(componentsFor(named1));

		final int[]  c2 = alpha(componentsFor(named2));
		
		return shadeOf(
				f,
				c1[0],c1[1],c1[2],c1[3],
				c2[0],c2[1],c2[2],c2[3]
				);
		}
	
	public static String shadeOf(
			float f,
			int ra, int ga, int ba,int aa,
			int rb, int gb, int bb,int ab
			)
		{
		return toRGB(
			(ra+(int)((rb-ra)*f)),
			(ga+(int)((gb-ga)*f)),
			(ba+(int)((bb-ba)*f)),
			(aa+(int)((ab-aa)*f))
			);
		}
	
	public static String shadeOf(
			float f,
			int ra, int ga, int ba,
			int rb, int gb, int bb
			)
		{
		return shadeOf(f,ra,ga,ba,255,rb,gb,bb,255);
		}
}
