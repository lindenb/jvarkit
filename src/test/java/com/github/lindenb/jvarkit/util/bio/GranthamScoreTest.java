package com.github.lindenb.jvarkit.util.bio;

import org.testng.Assert;
import org.testng.annotations.Test;

public class GranthamScoreTest {
@Test
public void test01() {
	Assert.assertEquals(GranthamScore.score('_', '_'), GranthamScore.getDefaultScore());
	final String nucl="ARNDCQEGHILKMFPSTWYV";
	for(int i=0;i+1< nucl.length();i++) {
		for(int j=i+1;j< nucl.length();j++) {
			if(i==j) 	Assert.assertEquals(GranthamScore.score(nucl.charAt(i),nucl.charAt(j)),0);
			
			
			Assert.assertEquals(
					GranthamScore.score(nucl.charAt(i),nucl.charAt(j)),
					GranthamScore.score(nucl.charAt(j),nucl.charAt(i))
					);
			
			Assert.assertTrue(GranthamScore.score(nucl.charAt(i),nucl.charAt(j))>=0);
			Assert.assertEquals(
					GranthamScore.score(nucl.charAt(i),nucl.charAt(j)),
					GranthamScore.score(Character.toLowerCase(nucl.charAt(i)),nucl.charAt(j))
					);
			Assert.assertEquals(
					GranthamScore.score(nucl.charAt(i),nucl.charAt(j)),
					GranthamScore.score(nucl.charAt(i),Character.toLowerCase(nucl.charAt(j)))
					);
			Assert.assertEquals(
					GranthamScore.score(nucl.charAt(i),nucl.charAt(j)),
					GranthamScore.score(Character.toLowerCase(nucl.charAt(i)),Character.toLowerCase(nucl.charAt(j)))
					);
		}
	}
	
	}
}
