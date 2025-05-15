package com.github.lindenb.jvarkit.variant.compact;

import java.nio.charset.Charset;
import java.util.Arrays;

class VCFCompactBase {
	protected static final Charset ENCODING= Charset.forName("UTF-8");
	interface LightGenotype {
		int[] getIndexes();
		default boolean isMissing() {
			int[] array=getIndexes();
			return array==null || array.length==0 || Arrays.stream(array).anyMatch(I->I==-1);
			}
		default boolean hasAlt() {
			int[] array=getIndexes();
			return array!=null || array.length>0 || Arrays.stream(array).anyMatch(I->I>0);
			}
		default int getPloidy() {
			int[] array=getIndexes();
			return array==null?0:array.length;
			}
		default boolean isHomRef() {
			int[] array=getIndexes();
			return array!=null && array.length==2 && array[0]==0 && array[1]==0;
			}
		default boolean isHomVar() {
			int[] array=getIndexes();
			return array!=null && array.length==2 && array[0]==1 && array[1]==1;
			}
		default boolean isHet() {
			int[] array=getIndexes();
			return array!=null && array.length==2 && (
					(array[0]==0 && array[1]==1) ||
					(array[0]==1 && array[1]==0)
					);
			}
		}
	public VCFCompactBase() {
		// TODO Auto-generated constructor stub
	}

}
