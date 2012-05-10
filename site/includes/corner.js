/** 
 * corner.js 1.82 (18-Mar-2010) (c) by Christian Effenberger 
 * All Rights Reserved. Source: corner.netzgesta.de
 * Distributed under Netzgestade Software License Agreement.
 * This license permits free of charge use on non-commercial 
 * and private web sites only under special conditions. 
 * Read more at... http://www.netzgesta.de/cvi/LICENSE.html	
 * Commercial licenses available via... cvi[at]netzgesta[dot]de	 
**/
//Modified by Richard Barnes to provide this service automagically to all iamges.
var tmp = navigator.appName == 'Microsoft Internet Explorer' && navigator.userAgent.indexOf('Opera') < 1 ? 1 : 0;
if(tmp) var isIE = document.namespaces ? 1 : 0;

if(isIE) {
	if(document.namespaces['v']==null) {
		var e=["shape","shapetype","group","background","path","formulas","handles","fill","stroke","shadow","textbox","textpath","imagedata","line","polyline","curve","roundrect","oval","rect","arc","image"],s=document.createStyleSheet(); 
		for(var i=0; i<e.length; i++) {s.addRule("v\\:"+e[i],"behavior: url(#default#VML);");} document.namespaces.add("v","urn:schemas-microsoft-com:vml");
	} 
}

function getRadius(radius,width,height){
	var part = (Math.min(width,height)/100);
	radius = Math.max(Math.min(100,radius/part),0);
	return radius + '%';
}

function getImages(className){
	var children = document.getElementsByTagName('img'); 
	var elements = new Array(); var i = 0;
	var child; var classNames; var j = 0;
	for (i=0;i<children.length;i++) {
		child = children[i];
		//classNames = child.className.split(' ');
		//for (var j = 0; j < classNames.length; j++) {
			//if (classNames[j] == className) {
				elements.push(child);
			//	break;
			//}
		//}
	}
	return elements;
}

function getClasses(classes,string){
	var temp = '';
	for (var j=0;j<classes.length;j++) {
		if (classes[j] != string) {
			if (temp) {
				temp += ' '
			}
			temp += classes[j];
		}
	}
	return temp;
}

function getClassValue(classes,string){
	var temp = 0; var pos = string.length;
	for (var j=0;j<classes.length;j++) {
		if (classes[j].indexOf(string) == 0) {
			temp = Math.min(classes[j].substring(pos),100);
			break;
		}
	}
	return Math.max(0,temp);
}

function getClassAttribute(classes,string){
	var temp = 0; var pos = string.length;
	for (var j=0;j<classes.length;j++) {
		if (classes[j].indexOf(string) == 0) {
			temp = 1; break;
		}
	}
	return temp;
}

function roundedRect(ctx,x,y,width,height,radius,nopath){
	if (!nopath) ctx.beginPath();
	ctx.moveTo(x,y+radius);
	ctx.lineTo(x,y+height-radius);
	ctx.quadraticCurveTo(x,y+height,x+radius,y+height);
	ctx.lineTo(x+width-radius,y+height);
	ctx.quadraticCurveTo(x+width,y+height,x+width,y+height-radius);
	ctx.lineTo(x+width,y+radius);
	ctx.quadraticCurveTo(x+width,y,x+width-radius,y);
	ctx.lineTo(x+radius,y);
	ctx.quadraticCurveTo(x,y,x,y+radius);
	if (!nopath) ctx.closePath();
}

function addGradient(ctx,x,y,w,h,color,opacity) {
	var tmp = ctx.createLinearGradient(x,y,w,h);
	var val = (color>0?0.25:0.2);
	tmp.addColorStop(0,'rgba('+color+','+color+','+color+',0.9)');
	tmp.addColorStop(val,'rgba('+color+','+color+','+color+','+opacity+')');
	tmp.addColorStop(0.75,'rgba('+color+','+color+','+color+',0)');
	tmp.addColorStop(1,'rgba('+color+','+color+','+color+',0)');
	return tmp;
}

function addRadialStyle(ctx,x1,y1,r1,x2,y2,r2,opacity) {
	var tmp = ctx.createRadialGradient(x1,y1,r1,x2,y2,r2);
	var opt = Math.min(parseFloat(opacity+0.1),1.0);
	tmp.addColorStop(0,'rgba(0,0,255,'+opt+')');
	tmp.addColorStop(0.25,'rgba(0,0,255,'+opacity+')');
	tmp.addColorStop(1,'rgba(0,0,255,0)');
	return tmp;
}

function addLinearStyle(ctx,x,y,w,h,opacity) {
	var tmp = ctx.createLinearGradient(x,y,w,h);
	var opt = Math.min(parseFloat(opacity+0.1),1.0);
	tmp.addColorStop(0,'rgba(0,0,255,'+opt+')');
	tmp.addColorStop(0.25,'rgba(0,0,255,'+opacity+')');
	tmp.addColorStop(1,'rgba(0,0,255,0)');
	return tmp;
}

function addShine(ctx,width,height,radius,opacity,extra) {
	var style; var color = (extra!=1?254:0);
	style = addGradient(ctx,0,radius,radius,radius,color,opacity);
	ctx.beginPath();
	ctx.moveTo(0,0);
	ctx.lineTo(0,height);	
	ctx.lineTo(radius,height);
	ctx.lineTo(radius,radius);
	ctx.closePath();
	ctx.fillStyle = style;
	ctx.fill();
	style = addGradient(ctx,radius,0,radius,radius,color,opacity);
	ctx.beginPath();
	ctx.moveTo(0,0);
	ctx.lineTo(width,0);
	ctx.lineTo(width,radius);
	ctx.lineTo(radius,radius);
	ctx.closePath();
	ctx.fillStyle = style;
	ctx.fill();
}

function addShade(ctx,width,height,radius,opacity) {
	var style;
	style = addGradient(ctx,width,radius,width-radius,radius,0,opacity);
	ctx.beginPath();
	ctx.moveTo(width,0);
	ctx.lineTo(width,height);
	ctx.lineTo(width-radius,height-radius);
	ctx.lineTo(width-radius,0);
	ctx.closePath();
	ctx.fillStyle = style;
	ctx.fill(); 
	style = addGradient(ctx,radius,height,radius,height-radius,0,opacity);
	ctx.beginPath();
	ctx.moveTo(width,height);
	ctx.lineTo(0,height);
	ctx.lineTo(0,height-radius);
	ctx.lineTo(width-radius,height-radius);
	ctx.closePath();
	ctx.fillStyle = style;
	ctx.fill();
}

function roundedShadow(ctx,x,y,width,height,radius,opacity){
	var style; 
	ctx.beginPath();
	ctx.rect(x,y+height-radius,radius,radius);
	ctx.closePath();
	style = addRadialStyle(ctx,x+radius,y+height-radius,radius-x,x+radius,y+height-radius,radius,opacity);
	ctx.fillStyle = style;
	ctx.fill();
	ctx.beginPath();
	ctx.rect(x+radius,y+height-y,width-(radius*2.25),y);
	ctx.closePath();
	style = addLinearStyle(ctx,x+radius,y+height-y,x+radius,y+height,opacity);
	ctx.fillStyle = style;
	ctx.fill();
	ctx.beginPath(); 
	ctx.rect(x+width-(radius*1.25),y+height-(radius*1.25),radius*1.25,radius*1.25);
	ctx.closePath();
	style = addRadialStyle(ctx,x+width-(radius*1.25),y+height-(radius*1.25),(radius*1.25)-1.5-x,x+width-(radius*1.25),y+height-(radius*1.25),radius*1.25,opacity);
	ctx.fillStyle = style;
	ctx.fill();
	ctx.beginPath();
	ctx.rect(x+width-x,y+radius,x,height-(radius*2.25));
	ctx.closePath();
	style = addLinearStyle(ctx,x+width-x,y+radius,x+width,y+radius,opacity);
	ctx.fillStyle = style;
	ctx.fill();
	ctx.beginPath();
	ctx.rect(x+width-radius,y,radius,radius);
	ctx.closePath();
	style = addRadialStyle(ctx,x+width-radius,y+radius,radius-x,x+width-radius,y+radius,radius,opacity);
	ctx.fillStyle = style;
	ctx.fill();
}

function addIECorners() {
	var theimages = getImages('corner');
	var image; var object; var vml;  var div; var pos; var i; var classes = '';
	var iradius = null; var ishadow = null; var ishade = null; var inverse = null; 
	var newClasses = ''; var maxdim = null; var offset = null; var radius = null; 
	var display = ""; var flt = null; var width = null; var height = null;
	var start, head, soft, shadow, fill, foot, end;
	var left, top, bottom, right, lt, br, linear, inset;
	for (i=0;i<theimages.length;i++) {	
		image = theimages[i];
		object = image.parentNode; 
		classes = image.className.split(' ');
		if (getClassAttribute(classes,"nobord")) { continue; }
		iradius = getClassValue(classes,"iradius");
		ishadow = getClassValue(classes,"ishadow");
		if (object.tagName=='A') { //getClassValue(classes,"ishade");
			ishadow = 50;
		} else {
			ishadow = 0;
		}
		ishade  = getClassValue(classes,"ishade");
		inverse = getClassAttribute(classes,"inverse");
		newClasses = getClasses(classes,"corner");
		width = image.width; height = image.height;
		maxdim = Math.min(width,height)/2;
		iradius = Math.min(maxdim,iradius); offset = 4;
		offset = (ishadow>0?(inverse>0?0:Math.min(Math.max(offset,iradius/2),16)):0);
		radius = getRadius(iradius,width,height);
		display = (image.currentStyle.display.toLowerCase()=='block')?'block':'inline-block';
		vml = document.createElement(['<var style="overflow:hidden;display:' + display + ';width:' + width + 'px;height:' + height + 'px;padding:0;">'].join(''));
		flt = image.currentStyle.styleFloat.toLowerCase();
		display = (flt=='left'||flt=='right')?'inline':display;
		start = '<v:group style="zoom:1; display:' + display + '; margin:-1px 0 0 -1px; padding:0; position:relative; width:' + width + 'px;height:' + height + 'px;" coordsize="' + width + ',' + height + '">';
		fill = '<v:fill src="' + image.src + '" type="frame" />';
		foot = (iradius>0?'</v:roundrect>':'</v:rect>'); end = '</v:group>';
		soft = ''; shadow = ''; lt = ''; left = ''; top = ''; bottom = ''; right = '';
		if(iradius<=0) {
			if (ishadow>0) {
				if(inverse<=0) {
					ishadow = ishadow/50; offset = 8; inset = 4;
					soft = '<v:rect strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0 ;width:' + width + 'px;height:' + height + 'px;"><v:fill color="#ffffff" opacity="0.0" /></v:rect><v:rect strokeweight="0" filled="t" stroked="f" fillcolor="#000000" style="filter:Alpha(opacity=' + (ishadow*64) + '), progid:dxImageTransform.Microsoft.Blur(PixelRadius=' + inset + ', MakeShadow=false); zoom:1;margin:-1px 0 0 -1px;padding: 0;display:block;position:absolute;top:' + inset + 'px;left:' + inset + 'px;width:' + (width-(3*inset)) + 'px;height:' + (height-(3*inset)) + 'px;"><v:fill color="#000000" opacity="1" /></v:rect>';
					head = '<v:rect strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0 ;width:' + (width-offset) + 'px;height:' + (height-offset) + 'px;">';
				}else if(inverse>0) {
					ishadow = ishadow/50; iradius = 12; linear = "linear";
					head = '<v:rect filled="t" stroked="t" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0; width:' + width + 'px;height:' + height + 'px;">';
					shadow = '<v:stroke weight="0.5" opacity="'+(ishadow/2)+'" color="#000000" />';
					top = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + width + ',' + iradius + '" path="m 0,0 l '+width+',0,'+width+','+iradius+','+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + width + 'px; height:' + iradius + 'px;"><v:fill method="' + linear + '" type="gradient" angle="0" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishadow + '" /></v:shape>'; 
					left = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + iradius + ',' + height + '" path="m 0,0 l 0,'+height+','+iradius+','+height+','+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + iradius + 'px; height:' + height + 'px;"><v:fill method="' + linear + '" type="gradient" angle="90" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishadow + '" /></v:shape>';
				}
			}else {
				head = '<v:rect strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="margin:-1px 0 0 -1px;padding:0;display:' + display + ';width:' + width + 'px;height:' + height + 'px;">';
			}
			if(ishade>0) {
				ishade = ishade/50; iradius = 12;
				top = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" coordorigin="0,0" coordsize="' + (width-offset) + ',' + iradius + '" path="m 0,0 l '+(width-offset)+',0,'+(width-offset)+','+iradius+','+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + (width-offset) + 'px; height:' + iradius + 'px;"><v:fill method="linear" type="gradient" angle="0" color="#ffffff" opacity="0" color2="#ffffff" o:opacity2="' + ishade + '" /></v:shape>'; 
				left = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" coordorigin="0,0" coordsize="' + iradius + ',' + (height-offset) + '" path="m 0,0 l 0,'+(height-offset)+','+iradius+','+(height-offset)+','+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + iradius + 'px; height:' + (height-offset) + 'px;"><v:fill method="linear" type="gradient" angle="90" color="#ffffff" opacity="0" color2="#ffffff" o:opacity2="' + ishade + '" /></v:shape>';
				bottom = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + (width-offset) + ',' + iradius + '" path="m 0,'+iradius+' l '+(width-offset)+','+iradius+','+(width-offset-iradius)+',0,'+iradius+',0 x e" style="position:absolute; margin: 0; top:' + (height-offset-iradius) + 'px; left: 0px; width:' + (width-offset) + 'px; height:' + iradius + 'px;"><v:fill method="linear" type="gradient" angle="180" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishade + '" /></v:shape>';
				right = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + iradius + ',' + (height-offset) + '" path="m '+iradius+',0 l '+iradius+','+(height-offset)+',0,'+(height-offset-iradius)+',0,'+iradius+' x e" style="position:absolute; margin: 0; top: 0px; left:' + (width-offset-iradius) + 'px; width:' + iradius + 'px; height:' + (height-offset) + 'px;"><v:fill method="linear" type="gradient" angle="270" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishade + '" /></v:shape>';
			}
		} else {
			if(ishadow>0) {
				linear = "linear sigma"; pos = 2;
				if(inverse<=0) {
					ishadow = ishadow/50; inset = Math.round(offset*0.5);
					soft = '<v:rect strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0 ;width:' + width + 'px;height:' + height + 'px;"><v:fill color="#ffffff" opacity="0.0" /></v:rect><v:roundrect arcsize="' + (radius+inset) + '" strokeweight="0" filled="t" stroked="f" fillcolor="#000000" style="filter:Alpha(opacity=' + (ishadow*64) + '), progid:dxImageTransform.Microsoft.Blur(PixelRadius=' + inset + ', MakeShadow=false); zoom:1;margin:-1px 0 0 -1px;padding: 0;display:block;position:absolute;top:' + inset + 'px;left:' + inset + 'px;width:' + (width-(3*inset)) + 'px;height:' + (height-(3*inset)) + 'px;"><v:fill color="#000000" opacity="1" /></v:roundrect>';
					head = '<v:roundrect arcsize="' + radius + '" strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0; width:' + (width-offset) + 'px;height:' + (height-offset) + 'px;">';
				}else if(inverse>0) {
					ishadow = ishadow/50;
					head = '<v:roundrect arcsize="' + radius + '" filled="t" stroked="t" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0; width:' + width + 'px;height:' + height + 'px;">';
					shadow = '<v:stroke weight="0.5" opacity="'+(ishadow/2)+'" color="#000000" />';
					top = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + width + ',' + iradius + '" path="m '+iradius+','+iradius+' l '+width+','+iradius+' qy '+(width-iradius)+',0 l '+iradius+',0 x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: -1px; width:' + (width+1) + 'px; height:' + iradius + 'px;"><v:fill method="' + linear + '" type="gradient" angle="0" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishadow + '" /></v:shape>'; 
					left = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + iradius + ',' + height + '" path="m 0,'+iradius+' l 0,'+(height-iradius)+' qy '+iradius+','+height+' l '+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: -1px; left: 0px; width:' + iradius + 'px; height:' + (height+1) + 'px;"><v:fill method="' + linear + '" type="gradient" angle="90" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishadow + '" /></v:shape>';
					lt = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#000000" coordorigin="0,0" coordsize="' + iradius + ',' + iradius + '" path="m '+iradius+','+iradius+' l 0,'+iradius+' qy '+iradius+',0 l '+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + iradius + 'px; height:' + iradius + 'px;"><v:fill method="' + linear + '" focus="1" focusposition="1,1" focussize="0.5,0.5" type="gradientradial" color="#000000" opacity="0" color2="#000000" o:opacity2="' + ishadow + '" /></v:shape>';
				}
			}else {
				pos = 1; offset = 0;
				head = '<v:roundrect arcsize="' + radius + '" strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" style="position:absolute; margin:-1px 0 0 -1px;padding:0; width:' + width + 'px;height:' + height + 'px;">';
			}
			if(ishade>0) {
				ishade = ishade/50; linear = "linear";
				top = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" coordorigin="0,0" coordsize="' + (width-offset) + ',' + iradius + '" path="m '+iradius+','+iradius+' l '+(width-offset)+','+iradius+' qy '+(width-offset-iradius)+',0 l '+iradius+',0 x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: -1px; width:' + (width-offset+pos) + 'px; height:' + iradius + 'px;"><v:fill method="' + linear + '" type="gradient" angle="0" color="#ffffff" opacity="0" color2="#ffffff" o:opacity2="' + ishade + '" /></v:shape>'; 
				left = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" coordorigin="0,0" coordsize="' + iradius + ',' + (height-offset) + '" path="m 0,'+iradius+' l 0,'+(height-iradius-offset)+' qy '+iradius+','+(height-offset)+' l '+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: -1px; left: 0px; width:' + iradius + 'px; height:' + (height-offset+pos) + 'px;"><v:fill method="' + linear + '" type="gradient" angle="90" color="#ffffff" opacity="0" color2="#ffffff" o:opacity2="' + ishade + '" /></v:shape>';
				lt = '<v:shape strokeweight="0" filled="t" stroked="f" fillcolor="#ffffff" coordorigin="0,0" coordsize="' + iradius + ',' + iradius + '" path="m '+iradius+','+iradius+' l 0,'+iradius+' qy '+iradius+',0 l '+iradius+','+iradius+' x e" style="position:absolute; margin: -1px 0 0 -1px; top: 0px; left: 0px; width:' + iradius + 'px; height:' + iradius + 'px;"><v:fill method="' + linear + '" focus="1" focusposition="1,1" focussize="0.5,0.5" type="gradientradial" color="#ffffff" opacity="0" color2="#ffffff" o:opacity2="' + ishade + '" /></v:shape>';
			}
		}
		vml.innerHTML = start + soft + head + fill + shadow + foot + right + bottom + top + left + lt + end;
		vml.className = newClasses;
		vml.style.cssText = image.style.cssText;
		vml.style.height = image.height+'px';
		vml.style.width = image.width+'px';
		vml.height = image.height;
		vml.width = image.width;
		vml.src = image.src; vml.alt = image.alt;
		if(image.id!='') vml.id = image.id; 
		if(image.title!='') vml.title = image.title;
		if(image.getAttribute('onclick')!='') vml.setAttribute('onclick',image.getAttribute('onclick'));
		if(image.getAttribute("usemap")) {
			if (iradius>0){pos = offset;}else {pos = 0;}
			object.style.position = 'relative';
			object.style.height = height+'px';
			object.style.width = width+'px';
			image.left = 0; image.top = 0;
			image.style.position = 'absolute';
			image.style.height = height+'px';
			image.style.width = width+'px';
			image.style.left = 0 + 'px';
			image.style.top = 0 + 'px';
			image.style.filter = "Alpha(opacity=0)";
			object.insertBefore(vml,image);
		}else {
			object.replaceChild(vml,image);
		}
	}
}

function addCorners() {
	var theimages = getImages('corner');
	var image; var object; var canvas; var context; var i;
	var iradius = null; var ishade = null; var ishadow = null;
	var inverse = null; var classes = ''; var newClasses = ''; 
	var maxdim = null; var style = null; var offset = null;
	for (i=0;i<theimages.length;i++) {
		image = theimages[i];
		object = image.parentNode; 
		canvas = document.createElement('canvas');
		if (canvas.getContext) {
			classes = image.className.split(' ');
			if (getClassAttribute(classes,"nobord")) { continue; }
			iradius = 15; //getClassValue(classes,"iradius");
			ishadow = getClassValue(classes,"ishadow");
			if (object.tagName=='A') { //getClassValue(classes,"ishade");
				ishadow = 50;
			} else {
				ishadow = 0;
			}
			inverse = getClassAttribute(classes,"inverse");
			newClasses = getClasses(classes,"corner");
			canvas.className = newClasses;
			canvas.style.cssText = image.style.cssText;
			canvas.style.height = image.height+'px';
			canvas.style.width = image.width+'px';
			canvas.height = image.height;
			canvas.width = image.width;
			canvas.src = image.src; canvas.alt = image.alt;
			if(image.id!='') canvas.id = image.id;
			if(image.title!='') canvas.title = image.title;
			if(image.getAttribute('onclick')!='') canvas.setAttribute('onclick',image.getAttribute('onclick'));
			maxdim = Math.min(canvas.width,canvas.height)/2;
			iradius = Math.min(maxdim,iradius); offset = 4;
			offset = (ishadow>0?(inverse>0?0:Math.min(Math.max(offset,iradius/2),16)):0);
			context = canvas.getContext("2d");
			if(image.getAttribute("usemap")) {
				object.style.position = 'relative';
				object.style.height = image.height+'px';
				object.style.width = image.width+'px';
				canvas.left = 0; canvas.top = 0;
				canvas.style.position = 'absolute';
				canvas.style.left = 0 + 'px';
				canvas.style.top = 0 + 'px';
				image.left = 0; image.top = 0;
				image.style.position = 'absolute';
				image.style.height = image.height+'px';
				image.style.width = image.width+'px';
				image.style.left = 0 + 'px';
				image.style.top = 0 + 'px';
				image.style.opacity = 0;
				object.insertBefore(canvas,image);
			}else {
				object.replaceChild(canvas,image);
			}
			context.clearRect(0,0,canvas.width,canvas.height);
			context.save();
			if (ishadow>0 && inverse<=0) {
				ishadow = ishadow/100;
				if (iradius>0) {
					roundedShadow(context,offset,offset,canvas.width-offset,canvas.height-offset,iradius,ishadow);
				}else {
					offset = 8; 
					roundedShadow(context,offset,offset,canvas.width-offset,canvas.height-offset,offset,ishadow);
				}
			}
			if (iradius<=0) {
				context.beginPath();
				context.rect(0,0,canvas.width-offset,canvas.height-offset);
				context.closePath();
			}else {
				roundedRect(context,0,0,canvas.width-offset,canvas.height-offset,iradius);
			}
			context.clip();
			context.fillStyle = 'rgba(0,0,0,0)';
			context.fillRect(0,0,canvas.width,canvas.height);
			context.drawImage(image,0,0,canvas.width-offset,canvas.height-offset);
			if (ishadow>0 && inverse>0) {
				ishadow = ishadow/100;
				if (iradius>0) {
					addShine(context,canvas.width,canvas.height,iradius,ishadow,1);
					roundedRect(context,0,0,canvas.width,canvas.height,iradius);
				}else {
					iradius = 16; 
					addShine(context,canvas.width,canvas.height,iradius,ishadow,1);
					context.beginPath();
					context.rect(0,0,canvas.width,canvas.height);
					context.closePath();
				}
				context.strokeStyle = 'rgba(0,0,0,'+ishadow+')';
				context.lineWidth = 2;
				context.stroke();
			}
			if (ishade>0) {
				ishade = ishade/100;
				if (iradius<=0) iradius = 16; 
				addShade(context,canvas.width-offset,canvas.height-offset,iradius,ishade);
				addShine(context,canvas.width-offset,canvas.height-offset,iradius,ishade);
			}
			canvas.style.visibility = 'visible';
		}
	}
}

if(window.attachEvent&&!window.opera) window.attachEvent("onload",addIECorners);
else window.addEventListener("load",addCorners,false);
