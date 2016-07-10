
var gbrowser = new GenomeBrowser();

	createCheckbox({"id":"properpair","text":"Proper Pairs","checked":true});
	createCheckbox({"id":"firstinpair","text":"First In Pair","checked":true});
	createCheckbox({"id":"secondinpair","text":"Second In Pair","checked":true});
	createCheckbox({"id":"","text":"Not primary alignment","checked":false});
	createCheckbox({"id":"","text":"Fails Quality Check","checked":false});
	createCheckbox({"id":"showdup","text":"Duplicates","checked":false});
	createCheckbox({"id":"","text":"Supplementary Align","checked":false});

function paintConfigIndex(idx) {
	if(idx<0 || idx >= config.length) return;

	var xmlhttp = new XMLHttpRequest();
	xmlhttp.onreadystatechange = function() {
		if (xmlhttp.readyState == 4 && xmlhttp.status == 200) {
		    var responseText = JSON.parse(xmlhttp.responseText);
		    var params = {"reads":[],"canvasid":"canvasdoc"};
		    var showdup = document.getElementById("showdup").checked;
		    var properpair = document.getElementById("properpair").checked;
		    var firstinpair = document.getElementById("firstinpair").checked;
		    var secondinpair = document.getElementById("secondinpair").checked;
		    var secondaryalign = document.getElementById("secondaryalign").checked;
		    var failsqc = document.getElementById("failsqc").checked;
		    var supalign = document.getElementById("supalign").checked;
		    
		      
		    
		    
			for(var i in responseText.reads) {
				var rec = new SamRecord(responseText.reads[i]);
				if( rec.isReadUnmappedFlag()) continue;
				/*
				if( !showdup && rec.getDuplicateReadFlag()) continue;
				if( !properpair && !rec.isProperPairFlag()) continue;
				if( firstinpair && !rec.getFirstInPairFlag()) continue;
				if( secondinpair && !rec.getSecondInPairFlag()) continue;
				if( failsqc && !rec.getReadFailsVendorQualityCheckFlag()) continue;
				if( supalign && rec.isSupplementaryAlignmentFlag()) continue;
				*/
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
	var e = document.createElement("label");
	e.setAttribute("for",cb.id);
	e.appendChild(document.createTextNode(cb.text));
	span.appendChild(e);
	e = document.createElement("input");
	span.appendChild(e);
	e.setAttribute("type","checkbox");
	e.setAttribute("id",cb.id);
	e.setAttribute("checked",cb.checked);
	e.addEventListener("change",repaintConfig,false);
	
}

function init()
	{
	if( document.getElementById("flags")==null) console.log("BOUUUMMMM");
	
	createCheckbox({"id":"properpair","text":"Only Proper Pairs","checked":false});
	createCheckbox({"id":"firstinpair","text":"Only First In Pair","checked":false});
	createCheckbox({"id":"secondinpair","text":"Only Second In Pair","checked":false});
	createCheckbox({"id":"secondaryalign","text":"Not primary alignment","checked":false});
	createCheckbox({"id":"failsqc","text":"Fails Quality Check","checked":false});
	createCheckbox({"id":"showdup","text":"Duplicates","checked":false});
	createCheckbox({"id":"supalign","text":"Supplementary Align","checked":false});
	
	
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
