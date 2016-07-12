
var gbrowser = new GenomeBrowser();


function paintConfigIndex(idx) {
	if(idx<0 || idx >= config.length) return;

	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
		    var responseText = JSON.parse(xmlhttp.responseText);
		    var params = {"reads":[],"canvasid":"canvasdoc",
		    	"interval": new Interval(responseText.interval),
				"header":null
		    	};
		    if( "reference" in  responseText)
		    	{
		    	params.reference = new ReferenceSequence(
		    		params.interval.contig,
		    		params.interval.start,
		    		responseText.reference
		    		);
		    	}
		    if( "title" in responseText)
		    	{
		    	document.getElementById("browserTitle").innerHTML = responseText.title;
		    	}
		    else
		    	{
		    	document.getElementById("browserTitle").innerHTML = params.interval.toString();
		    	}
		    var HIDE=1;
		    var ONLY=2;
		    var duplicate = document.getElementById("duplicate").selectedIndex;
		    var shownonproperpair = document.getElementById("shownonproperpair").selectedIndex;
		    var firstinpair = document.getElementById("firstinpair").selectedIndex;
		    var secondinpair = document.getElementById("secondinpair").selectedIndex;
		    var secondaryalign = document.getElementById("secondaryalign").selectedIndex;
		    var failsqc = document.getElementById("failsqc").selectedIndex;
		    var supalign = document.getElementById("supalign").selectedIndex;
		    var minmapq = parseInt(document.getElementById("mapq").value);
		    var plusstrand = document.getElementById("plusstrand").selectedIndex;
		    var minusstrand = document.getElementById("minusstrand").selectedIndex;
		    var cigarregex = document.getElementById("cigarregex").value.trim(); cigarregex=(cigarregex==""?null:new RegExp(cigarregex));
		    var discordantchr = document.getElementById("discordantchr").selectedIndex;
		    gbrowser.useClip =  document.getElementById("showclip").checked;
		    gbrowser.expandInsertions =  document.getElementById("expandinsert").checked;
		    gbrowser.printReadBases =  document.getElementById("showreadbases").checked;
		    gbrowser.printReadName =  document.getElementById("printreadname").checked;
		    gbrowser.expandDeletions =  document.getElementById("expandeletion").checked;
		    
		    var responsereads;
		    if( "reads" in responseText)
				{
				responsereads = responseText.reads;
				}
		    else if( "sam" in responseText) {
				params.header = new SAMFileHeader(responseText.sam.header);
				responsereads = responseText.sam.reads;
				}
		    else {
				responsereads=[];
				}


			for(var i in responsereads) {
				var rec = new SamRecord(responsereads[i]);
				if( rec.isReadUnmappedFlag()) continue;
				
				if( rec.getMappingQuality() < minmapq)  continue;
				
				if(plusstrand>0) {
					if( plusstrand==HIDE && !rec.isReadPositiveStrandFlag()) continue;
					if( plusstrand==ONLY && rec.isReadPositiveStrandFlag()) continue;
				}
				
				if(minusstrand>0) {
					if( minusstrand==HIDE && !rec.isReadNegativeStrandFlag()) continue;
					if( minusstrand==ONLY && rec.isReadNegativeStrandFlag()) continue;
					}
				
				if( duplicate > 0 )
					{
					if( duplicate==HIDE && rec.getDuplicateReadFlag()) continue;
					if( duplicate==ONLY && !rec.getDuplicateReadFlag()) continue;
					}
				
				if( shownonproperpair > 0)
					{
					if( shownonproperpair==HIDE && !rec.isProperPairFlag()) continue;
					if( shownonproperpair==ONLY &&  rec.isProperPairFlag()) continue;
					}
				
				if( secondaryalign > 0)
					{
					if( secondaryalign==HIDE && rec.isNotPrimaryAlignmentFlag()) continue;
					if( secondaryalign==ONLY && !rec.isNotPrimaryAlignmentFlag()) continue;
					}
				
				if( failsqc > 0)
					{
					if( failsqc ==HIDE && rec.getReadFailsVendorQualityCheckFlag()) continue;
					if( failsqc ==ONLY && !rec.getReadFailsVendorQualityCheckFlag()) continue;
					}
				
				if( firstinpair > 0)
					{
					if( firstinpair ==HIDE && rec.getFirstInPairFlag()) continue;
					if( firstinpair ==ONLY && !rec.getFirstInPairFlag()) continue;
					}
				
				if( secondinpair > 0)
					{
					if( secondinpair ==HIDE && rec.getSecondInPairFlag()) continue;
					if( secondinpair ==ONLY && !rec.getSecondInPairFlag()) continue;
					}
				if( supalign > 0 )
					{
					if( supalign ==HIDE && rec.isSupplementaryAlignmentFlag()) continue;
					if( supalign ==ONLY && !rec.isSupplementaryAlignmentFlag()) continue;
					}
				if( cigarregex !=null && rec.getCigarString().search(cigarregex)==-1)
					{
					continue;
					}
					
				if(discordantchr>0)
					{
					var is_discordant= rec.hasDiscordantContigs();
					if( is_discordant && discordantchr==HIDE) continue;
					if( !is_discordant && discordantchr==ONLY) continue;
					}
				
				params.reads.push(rec);
				}
				
			try { gbrowser.paint(params); } catch(e) {alert(e+" "+e.stack);}
			}
		};
	xmlhttp.open("GET", config[idx].href, true);
	/** https://developer.mozilla.org/en-US/docs/Web/API/XMLHttpRequest/Using_XMLHttpRequest */
	xmlhttp.overrideMimeType("application/json;charset=UTF-8");
	xmlhttp.send();
	}

function repaintConfig() {
	var sel= document.getElementById("menu");
	paintConfigIndex(sel.selectedIndex);
	}

function changemenu(shift) {
	var sel= document.getElementById("menu");
	var s = sel.selectedIndex;
	s+=shift;
	if( s>= config.length) s=0;
	if( s<0) s = config.length-1;
	sel.selectedIndex = s;
	repaintConfig();
	}

function createCheckbox(cb) {
	var div = document.getElementById("flags");
	if( div == null) { console.log(cb.text+" "+cb.id); return; }
	var span = document.createElement("span");
	div.appendChild(span);
	
	var e = document.createElement("input");
	span.appendChild(e);
	e.setAttribute("type","checkbox");
	e.setAttribute("title",cb.text);
	e.setAttribute("id",cb.id);
	if(cb.checked) e.setAttribute("checked","true");
	e.setAttribute("value",cb.checked);
	e.addEventListener("change",repaintConfig,false);
	
	e = document.createElement("label");
	e.setAttribute("for",cb.id);
	e.appendChild(document.createTextNode(cb.text+". "));
	span.appendChild(e);
}

function createTextField(cb) {
	var div = document.getElementById("flags");
	var span = document.createElement("span");
	div.appendChild(span);
	var e = document.createElement("label");
	e.setAttribute("for",cb.id);
	e.appendChild(document.createTextNode(cb.text+":"));
	span.appendChild(e);
	e = document.createElement("input");
	span.appendChild(e);
	e.setAttribute("type","text");
	e.setAttribute("id",cb.id);
	
	e.setAttribute("value",cb.value);
	e.addEventListener("change",repaintConfig,false);
}


function createShowHide(cb) {
	var e,x,opts=["*","Hide","Only"]
	var div = document.getElementById("flags");
	var span = document.createElement("span");
	div.appendChild(span);
	
	e = document.createElement("select");
	span.appendChild(e);
	e.setAttribute("title",cb.text);
	e.setAttribute("id",cb.id);
	for(x in opts)
		{
		var o = document.createElement("option");
		o.setAttribute("value",opts[x]);
		if( (("index" in cb) && cb.index == x ) || (!("index" in cb) && x==0)) {
			o.setAttribute("selected","true");
			}
		o.appendChild(document.createTextNode(opts[x]+" "+cb.text));
		e.appendChild(o);
		}
	e.addEventListener("change",repaintConfig,false);
	}

function init()
	{
	if( document.getElementById("flags")==null) console.log("BOUUUMMMM");
	
	createShowHide({"id":"shownonproperpair","text":"Non Proper-Pairs","index":1});
	createShowHide({"id":"firstinpair","text":"First In Pair","index":0});
	createShowHide({"id":"secondinpair","text":"Second In Pair","index":0});
	createShowHide({"id":"secondaryalign","text":"Not primary alignment","index":1});
	createShowHide({"id":"failsqc","text":"Fails Quality Check","index":1});
	createShowHide({"id":"duplicate","text":"Duplicates","index":1});
	createShowHide({"id":"supalign","text":"Supplementary Align","index":1});
	createShowHide({"id":"plusstrand","text":"Strand (+)","index":0});
	createShowHide({"id":"minusstrand","text":"Strand (-)","index":0});
	//
	createCheckbox({"id":"showclip","text":"Show Clipped Regions","checked":false});
	createCheckbox({"id":"expandinsert","text":"Expand Insertion","checked":false});
	createCheckbox({"id":"expandeletion","text":"Expand Deletions","checked":true});
	
	createCheckbox({"id":"showreadbases","text":"Show Read Bases","checked":true});
	createTextField({"id":"cigarregex","text":"Cigar Regex","value":""});
	createCheckbox({"id":"printreadname","text":"Print Read Name","checked":false});
	createShowHide({"id":"discordantchr","text":"Discordant Contigs","index":1});
	
	
	createTextField({"id":"mapq","text":"Min Mapq","value":0});
	
	
	var menu= document.getElementById("menu");
	menu.addEventListener("change",function() {
		paintConfigIndex(menu.selectedIndex);
		},false);
	for( var i in config)
		{
		var cgf =  config[i];
		var opt = document.createElement("option");
		opt.setAttribute("value",i);
		menu.appendChild(opt);
		opt.appendChild(document.createTextNode(cgf.title));
		}
	paintConfigIndex(0);
	}


window.addEventListener("load",init,false);
