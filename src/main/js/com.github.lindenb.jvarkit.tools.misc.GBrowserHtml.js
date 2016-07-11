
var gbrowser = new GenomeBrowser();


function paintConfigIndex(idx) {
	if(idx<0 || idx >= config.length) return;

	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
		    var responseText = JSON.parse(xmlhttp.responseText);
		    var params = {"reads":[],"canvasid":"canvasdoc",
		    	"interval": new Interval(responseText.interval)
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
		    var showdup = document.getElementById("showdup").checked;
		    var shownonproperpair = document.getElementById("shownonproperpair").checked;
		    var hidefirstinpair = document.getElementById("hidefirstinpair").checked;
		    var hidesecondinpair = document.getElementById("hidesecondinpair").checked;
		    var secondaryalign = document.getElementById("secondaryalign").checked;
		    var showfailsqc = document.getElementById("showfailsqc").checked;
		    var showsupalign = document.getElementById("showsupalign").checked;
		    var minmapq = parseInt(document.getElementById("mapq").value);
		    var showplusstrand = document.getElementById("showplusstrand").checked;
		    var showminusstrand = document.getElementById("showminusstrand").checked;
		    var hidcigarpurem = document.getElementById("hidcigarpurem").checked;
		    var discordantchr = document.getElementById("discordantchr").selectedIndex;
		    gbrowser.useClip =  document.getElementById("showclip").checked;
		    gbrowser.expandInsertions =  document.getElementById("expandinsert").checked;
		    gbrowser.printReadBases =  document.getElementById("showreadbases").checked;
		    gbrowser.printReadName =  document.getElementById("printreadname").checked;
		    gbrowser.expandDeletions =  document.getElementById("expandeletion").checked;
		    
			for(var i in responseText.reads) {
				var rec = new SamRecord(responseText.reads[i]);
				if( rec.isReadUnmappedFlag()) continue;
				
				if( rec.getMappingQuality() < minmapq)  continue;
				if( !showplusstrand && rec.isReadPositiveStrandFlag()) continue;
				if( !showminusstrand && rec.isReadNegativeStrandFlag()) continue;
				
				if( !showdup && rec.getDuplicateReadFlag()) continue;
				
				if( !shownonproperpair && !rec.isProperPairFlag()) continue;
				
				if( hidefirstinpair && rec.getFirstInPairFlag()) continue;
				if( hidesecondinpair && rec.getSecondInPairFlag()) continue;
				if( !showfailsqc && rec.getReadFailsVendorQualityCheckFlag()) continue;
				if( !showsupalign && rec.isSupplementaryAlignmentFlag()) continue;
				if( hidcigarpurem && rec.getCigar().getNumElements()==1  && rec.getCigar().get(0).getOperator().isOneOf("M=")) continue;
				if(discordantchr>0)
					{
					var is_discordant= rec.hasDiscordantContigs();
					if( is_discordant && discordantchr==1) continue;
					if( !is_discordant && discordantchr==2) continue;
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
	var e,x,opts=["","Hide","Only"]
	var div = document.getElementById("flags");
	var span = document.createElement("span");
	div.appendChild(span);
	
	e = document.createElement("label");
	e.setAttribute("for",cb.id);
	e.appendChild(document.createTextNode(cb.text+":"));
	span.appendChild(e);
	
	e = document.createElement("select");
	span.appendChild(e);
	e.setAttribute("title",cb.text);
	e.setAttribute("id",cb.id);
	for(x in opts)
		{
		var o = document.createElement("option");
		e.setAttribute("value",opts[x]);
		if(x==0) e.setAttribute("selected","true");
		o.appendChild(document.createTextNode(opts[x]));
		e.appendChild(o);
		}
	e.addEventListener("change",repaintConfig,false);
}

function init()
	{
	if( document.getElementById("flags")==null) console.log("BOUUUMMMM");
	
	createCheckbox({"id":"shownonproperpair","text":"Show non Proper-Pairs","checked":false});
	createCheckbox({"id":"hidefirstinpair","text":"Hide First In Pair","checked":false});
	createCheckbox({"id":"hidesecondinpair","text":"Hide Second In Pair","checked":false});
	createCheckbox({"id":"secondaryalign","text":"Not primary alignment","checked":false});
	createCheckbox({"id":"showfailsqc","text":"Show Fails Quality Check","checked":false});
	createCheckbox({"id":"showdup","text":"Show Duplicates","checked":false});
	createCheckbox({"id":"showsupalign","text":"Show Supplementary Align","checked":false});
	createCheckbox({"id":"showplusstrand","text":"Show Strand +","checked":true});
	createCheckbox({"id":"showminusstrand","text":"Show Strand -","checked":true});
	//
	createCheckbox({"id":"showclip","text":"Show Clipped Regions","checked":false});
	createCheckbox({"id":"expandinsert","text":"Expand Insertion","checked":false});
	createCheckbox({"id":"expandeletion","text":"Expand Deletions","checked":true});
	
	createCheckbox({"id":"showreadbases","text":"Show Read Bases","checked":true});
	createCheckbox({"id":"hidcigarpurem","text":"Hide pure align","checked":false});
	createCheckbox({"id":"printreadname","text":"Print Read Name","checked":false});
	createShowHide({"id":"discordantchr","text":"Discordant Contigs"});
	
	
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
