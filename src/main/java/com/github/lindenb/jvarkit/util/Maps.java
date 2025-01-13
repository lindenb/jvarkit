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
package com.github.lindenb.jvarkit.util;

import java.util.LinkedHashMap;
import java.util.Map;


public class Maps {
	public static <K,V>  Map<K,V> of(K k1, V v1) {
		final Map<K,V> m=new LinkedHashMap<>();
		m.put(k1,v1);
		return m;
		}
	
	/** 
	 awk '{N=int($1);printf("public static <K,V>  Map<K,V> of(K k1, V v1");for(i=2;i<=N;i++) printf(",K k%d,V v%d",i,i); printf(") {\n\tfinal Map<K,V> m= of("); for(i=1;i<N;i++) printf("%sk%d,v%d",(i>1?",":""),i,i); printf("); m.put(k%d,v%d); return m;\n\t}\n",$1,$1);}'
	 */

	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2) {
		final Map<K,V> m= of(k1,v1); m.put(k2,v2); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3) {
		final Map<K,V> m= of(k1,v1,k2,v2); m.put(k3,v3); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3); m.put(k4,v4); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4); m.put(k5,v5); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5,K k6,V v6) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4,k5,v5); m.put(k6,v6); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5,K k6,V v6,K k7,V v7) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4,k5,v5,k6,v6); m.put(k7,v7); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5,K k6,V v6,K k7,V v7,K k8,V v8) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4,k5,v5,k6,v6,k7,v7); m.put(k8,v8); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5,K k6,V v6,K k7,V v7,K k8,V v8,K k9,V v9) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4,k5,v5,k6,v6,k7,v7,k8,v8); m.put(k9,v9); return m;
		}
	public static <K,V>  Map<K,V> of(K k1, V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5,K k6,V v6,K k7,V v7,K k8,V v8,K k9,V v9,K k10,V v10) {
		final Map<K,V> m= of(k1,v1,k2,v2,k3,v3,k4,v4,k5,v5,k6,v6,k7,v7,k8,v8,k9,v9); m.put(k10,v10); return m;
		}
	}
