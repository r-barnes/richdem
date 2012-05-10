<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
                version="1.0">

<xsl:output method="html"/>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!-- MdDlgHelp.xsl                                          -->
<!--                                                        -->    
<!-- Style sheet used to generate the HTML file containing  --> 
<!-- the function/variable help information.                -->
<!--                                                        -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                       Variables                        -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->

<xsl:variable name="pre_path" select="//MdElementDialogInfo/HelpPath/text()"/>
<xsl:variable name="path" select="translate($pre_path, '\', '/')"/>

<xsl:variable name="Font">arial,verdana</xsl:variable>
<xsl:variable name="Margin">0em 0em 0em 0em</xsl:variable>
<xsl:variable name="BackgroundColor">White</xsl:variable>
<xsl:variable name="SmallFont">verdana,arial</xsl:variable>
<xsl:variable name="SmallSize">7pt</xsl:variable>
<xsl:variable name="SmallColor">Black</xsl:variable>
<xsl:variable name="SmallWeight">Normal</xsl:variable>
<xsl:variable name="SmallStyle">Normal</xsl:variable>
<xsl:variable name="ParaFont">arial,verdana</xsl:variable>
<xsl:variable name="ParaSize">10pt</xsl:variable>
<xsl:variable name="ParaColor">Black</xsl:variable>
<xsl:variable name="ParaWeight">Normal</xsl:variable>
<xsl:variable name="ParaStyle">Normal</xsl:variable>
<xsl:variable name="ULFont">arial,verdana</xsl:variable>
<xsl:variable name="ULSize">10pt</xsl:variable>
<xsl:variable name="ULColor">Black</xsl:variable>
<xsl:variable name="ULWeight">Normal</xsl:variable>
<xsl:variable name="ULStyle">Normal</xsl:variable>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                 <ElementHelp> Template                 -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="ElementHelp">
  <p><xsl:apply-templates select="text()" /></p>
  <xsl:apply-templates select="*" />
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                    <para> Template                     -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="para">
	<p>
		<xsl:apply-templates select="*|text()"/>
	</p>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--               <bullet_item> Template                   -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="bullet_item">
	<ul>
		<li>
			<xsl:apply-templates select="*|text()"/>
		</li>
	</ul>
</xsl:template>


<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                  <indent> Template                     -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="indent">
	<dl>
		<dt/>
		<dd>
			<xsl:apply-templates select="*" />			
		</dd>
	</dl>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                    <link> Template                     -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="link">
	<a TARGET="viewer">
		<xsl:attribute name="href">
			<xsl:choose>
				<xsl:when test="contains(@src, 'ARCTOOLBOXHELP')">
						<xsl:call-template name="ToolPath">
							<xsl:with-param name="fullpath" select="@src"/>
						</xsl:call-template> 
				</xsl:when>
					<xsl:otherwise>
						<xsl:value-of select="@src"/>
					</xsl:otherwise>						
				</xsl:choose>
		</xsl:attribute>
		<xsl:value-of select="." />
	</a>		
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                   <illust> Template                    -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="illust">
<br/>
	<img>
		<xsl:attribute name="src">
			<xsl:choose>
				<xsl:when test="contains(@src, 'ARCTOOLBOXHELP')">
						<xsl:call-template name="ToolPath">
							<xsl:with-param name="fullpath" select="@src"/>
						</xsl:call-template> 
				</xsl:when>
				<xsl:otherwise>
						<xsl:value-of select="@src"/>
				</xsl:otherwise>						
			</xsl:choose>
		</xsl:attribute>
		
		<xsl:attribute name="alt">
			<xsl:value-of select="@alt"/>
		</xsl:attribute>
	</img>
	<p/>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                    <code> Template                     -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="code">
	<xsl:call-template name="codeMacro"/>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                 <codeMacro> Template                   -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template name ="codeMacro">
	<pre>
		<xsl:apply-templates select="*|text()"/>
	</pre>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                 <subSection> Template                  -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="subSection">
  <xsl:if test="*">
	  <p>
		  <h2>
			  <xsl:call-template name="ImageAndScript"/>
				  <xsl:text> </xsl:text><xsl:value-of select="@title"/>
		  </h2>	
	  </p>
	  <div>
		  <xsl:attribute name="Class">expand</xsl:attribute>
		  <xsl:attribute name="id"><xsl:value-of select="generate-id(.)"/></xsl:attribute>
		  <xsl:attribute name="style">display:None</xsl:attribute>
		  <xsl:apply-templates select="*"/>
	  </div>
  </xsl:if>
</xsl:template>	

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                <ImageAndScript> Template               -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template name="ImageAndScript">
  <a>
	  <xsl:variable name="NodeID" select="generate-id(.)"/>
		  <xsl:attribute name="href">#<xsl:value-of select="$NodeID"/></xsl:attribute>
		  <xsl:attribute name="onClick">expandIt('<xsl:value-of select="$NodeID"/>','<xsl:choose>
			  <xsl:when test="$path !=''">
				  <xsl:value-of select="concat($path, '/')"></xsl:value-of>	
			  </xsl:when>
			  <xsl:otherwise>
				  <xsl:value-of select="$path"></xsl:value-of>	
			  </xsl:otherwise>
		  </xsl:choose><xsl:text>','small_arrow_up.gif')</xsl:text>
		  </xsl:attribute>
		  <img width="11" height="11" alt="expand/collapse item" border="0" name="imEx" id="toolbox">
			  <xsl:attribute name="src">
					  <xsl:value-of select="$path"/>
						  <xsl:if test="$path !=''">
							  <xsl:text>/</xsl:text>
						  </xsl:if>
				  <xsl:text>small_arrow_up.gif</xsl:text>
			  </xsl:attribute>
		  </img>
  </a>
</xsl:template>	

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                   <ToolPath> Template                  -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template name="ToolPath">
	<xsl:param name="fullpath"/>
		<xsl:variable name="afterArcToolBoxHelp">
      <xsl:variable name="tempPath" select="substring-after($fullpath, '/')"/>
      <xsl:choose>
	      <xsl:when test="$tempPath !=''">
          <xsl:value-of select="$tempPath"/>
			  </xsl:when>
        <xsl:otherwise>
				  <xsl:value-of select="substring-after($fullpath, '\')"/>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:variable>

		<xsl:variable name="srcPath">	
			<xsl:value-of select="$path"/>
			<xsl:if test="$path !=''">
				<xsl:text>/</xsl:text>
			</xsl:if>
			<xsl:value-of select="$afterArcToolBoxHelp"/>
		</xsl:variable>
		<xsl:value-of select="$srcPath"></xsl:value-of>
</xsl:template>		

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--              MdElementDialogInfo Template              -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="MdElementDialogInfo">
<HTML>

<HEAD>
<TITLE>Help</TITLE>

<STYLE TYPE="text/css">

BODY {    font-family:      <xsl:value-of select="$Font"/>; 
          margin:           <xsl:value-of select="$Margin"/>; 
          background-color: <xsl:value-of select="$BackgroundColor"/>; }

.small {  font-family:      <xsl:value-of select="$SmallFont"/>; 
          font-size:        <xsl:value-of select="$SmallSize"/>; 
          color:            <xsl:value-of select="$SmallColor"/>; 
          font-weight:      <xsl:value-of select="$SmallWeight"/>; 
          font-style:       <xsl:value-of select="$SmallStyle"/>; }  

P {       font-family:      <xsl:value-of select="$ParaFont"/>; 
          font-size:        <xsl:value-of select="$ParaSize"/>; 
          color:            <xsl:value-of select="$ParaColor"/>; 
          font-weight:      <xsl:value-of select="$ParaWeight"/>; 
          font-style:       <xsl:value-of select="$ParaStyle"/>; }

UL {      font-family:      <xsl:value-of select="$ULFont"/>; 
          font-size:        <xsl:value-of select="$ULSize"/>; 
          color:            <xsl:value-of select="$ULColor"/>; 
          font-weight:      <xsl:value-of select="$ULWeight"/>; 
          font-style:       <xsl:value-of select="$ULStyle"/>; }

</STYLE><xsl:text>&#10;</xsl:text>

<xsl:comment> ================ Scripts ================ </xsl:comment>
<xsl:text>&#10;</xsl:text>

<SCRIPT src="{CommonPath}/depress.js" language="JavaScript"></SCRIPT><xsl:text>&#10;</xsl:text>

<SCRIPT language="JavaScript"><xsl:text>&#10;</xsl:text>
<xsl:comment>

<![CDATA[
function ShowArcGISHelp()
{
  window.external.ShowHelp();
}

function ShowHelpTopic(topic)
{
  document.getElementById("Intro").style.display = "none";
    
]]>
<xsl:for-each select="Properties/PropertyGroup"><xsl:for-each select="Property">  document.getElementById(&quot;<xsl:value-of select="PropertyName"/>Topic&quot;).style.display = "none";
</xsl:for-each></xsl:for-each><![CDATA[

  if (topic== '')
    topic = 'Intro';

  document.getElementById(topic).style.display = "block";
}

]]>
</xsl:comment><xsl:text>&#10;</xsl:text>
</SCRIPT><xsl:text>&#10;</xsl:text>

</HEAD><xsl:text>&#10;</xsl:text>

<BODY><xsl:text>&#10;</xsl:text>
<SCRIPT language="Javascript" src="{CommonPath}/disableclick.js"></SCRIPT><xsl:text>&#10;</xsl:text>

<xsl:comment>View Documentation...</xsl:comment><xsl:text>&#10;</xsl:text>

<xsl:text>&#10;</xsl:text>
<xsl:comment>General Description...</xsl:comment><xsl:text>&#10;</xsl:text>
<xsl:text>&#10;</xsl:text>

<DIV STYLE="margin:0.5em 0.5em 0.5em 1.0em;" ID="Intro" style="display:block;">
<B><xsl:value-of select="Title"/></B>
<xsl:apply-templates select="ElementHelp" />
</DIV><xsl:text>&#10;</xsl:text>

<xsl:text>&#10;</xsl:text>
<xsl:comment>Property Descriptions...</xsl:comment><xsl:text>&#10;</xsl:text>
<xsl:text>&#10;</xsl:text>

<xsl:for-each select="Properties/PropertyGroup">
  <xsl:for-each select="Property">
    <DIV STYLE="margin:0.5em 0.5em 0.5em 1.0em;" ID="{PropertyName}Topic" style="display:none;">
      <B><xsl:value-of select="PropertyLabel"/></B>
			<xsl:apply-templates select="PropertyHelp"/>
    </DIV><xsl:text>&#10;</xsl:text>
  </xsl:for-each>
</xsl:for-each>

<xsl:text>&#10;</xsl:text>
</BODY>
</HTML>
</xsl:template>

<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<!--                   Root Template                        -->
<!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ -->
<xsl:template match="/">
  <xsl:apply-templates />
</xsl:template>

</xsl:stylesheet>


