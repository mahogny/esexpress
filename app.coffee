#datasets = ["2i", "a2i", "lif", "s_lblast", "s_mblast", "s_eblast"]

GENESETTAB_CELLCOR      = "cellcor"
GENESETTAB_BETWEENGENES = "betweengenes"

tabcolorInactive = "#0033CC"
tabcolorActive   = "#9933CC"

publication_url = "http://example.com/new_url"

query = 
  gene: ""
dispidcounter=1

thecol = {}
thecol["mES_serum"] = "#AA0000"
thecol["mES_a2i"]   = "#EEAD0E"
thecol["mES_2i"]    = "#00008B"

thecol["sandberg_earlyblast"] = "#9ACD32"
thecol["sandberg_midblast"]   = "#789C31"  
thecol["sandberg_lateblast"]  = "#556B2F"

thecol["mES_2i: blastocyst-like"]    = "#0099CC"
thecol["mES_2i: 2C-like"]            = "#003399"
thecol["mES_serum: more pluripotent cells"]  = "#990000"
thecol["mES_serum: primed cells"]            = "#CC4D4D"
thecol["mES_serum: differentiating cells"]   = "#FF9999"



defaultshowds = []
defaultshowds["mES_serum"] = 1
defaultshowds["mES_a2i"] = 1
defaultshowds["mES_2i"]  = 1

mouseX=0
mouseY=0

current_div=null
current_divname=null

log2 = (val) ->
  Math.log(val) / Math.LN2;

###########################################################################
# Searching for genes
###########################################################################

search_genes_form = (e) ->
  e.preventDefault()
  search_genes()

search_genes = () ->
  query_url = "searchgene.php?gene=#{query.gene}"
  req = $.getJSON query_url
  req.success (data) ->
    search_gene_disp(data)
  req.fail (data) -> 
    alert("failed to query data, "+query_url)



search_gene_disp = (data) ->
  table = $ "#searchresult"
  table.empty()
  for rec in data
    link = $ "<a href=\"#\">"
    link.append rec.genesym
    link.click (() ->
      thisrec=rec
      return () ->
        if current_divname=="geneexp"
          (makeviewgenefunc thisrec.geneid, true)()
        else
          (makeviewgenefunc thisrec.geneid, false)())()
    table.append link
    table.append ($ "<br/>")



makegenelink = (geneid, genesym, showexplevel) ->
  link = $ "<a href=\"#\">"
  link.append genesym
  link.click (makeviewgenefunc geneid, showexplevel)

makeviewgenefunc = (geneid, showexplevel) ->
  return () ->
    hide_startpage()
    closecurrentdiv()
    view_gene(geneid, showexplevel)






###########################################################################
# Open up a single gene view
###########################################################################

view_gene = (gene, showexplevel) ->
  query.gene = gene
  query_url = "geneinfo.php?q=#{JSON.stringify query}"   
  req = $.getJSON query_url
  req.success (data) ->
    view_gene_disp(data, showexplevel)
  req.fail (data) -> 
    alert("failed to query data, "+query_url)


addtogeneset = (genesym) ->
  genesetlist = $ "#genesetlist"
  genesetlist.val (genesetlist.val()+"\n"+genesym).trim()


view_gene_disp = (geneinfo, showexplevel) ->
  dispidcounter++
  thisid="genedisp"+dispidcounter

  root = $ "#maincontent"
  form2 = $ "#genepanel"
  form2 = form2.html()
  form2 = form2.replace("GENEID", geneinfo.geneid)
  form2 = form2.replace("GENESYM", geneinfo.genesym)
  root.prepend form2
  form2 = root.children(":first")
  form2.attr id: thisid

  curtab = "geneexp"
  if !showexplevel
    curtab="corr"

  ##### Buttons for tabs
  updatetab = () ->
    form2.find("#tab_geneexp").addClass "hideclass"
    form2.find("#tab_genecorr").addClass "hideclass"
    form2.find("#btab_geneexp").css  "background-color", tabcolorInactive
    form2.find("#btab_genecorr").css "background-color", tabcolorInactive
    if(curtab=="geneexp")
      form2.find("#btab_geneexp").css "background-color", tabcolorActive
      form2.find("#tab_geneexp").removeClass "hideclass"
    else
      form2.find("#btab_genecorr").css "background-color", tabcolorActive
      form2.find("#tab_genecorr").removeClass "hideclass"
  updatetab()

  form2.find("#btab_genecorr").click () ->
    curtab="genecorr"
    updatetab()
  form2.find("#btab_geneexp").click () ->
    curtab="geneexp"
    updatetab()

  $("html, body").animate (scrollTop: form2.offset().top), "slow"

  form2.find("#genepanel-close").click () ->
    form2.remove()

  form2.find("#linkpubmed").attr href: ("http://www.ncbi.nlm.nih.gov/pubmed/?term="+geneinfo.genesym)
  form2.find("#linkensembl").attr href: ("http://www.ensembl.org/Mus_musculus/Gene/Summary?g="+geneinfo.geneid)
  form2.find("#linkgenecards").attr href: ("http://www.genecards.org/cgi-bin/carddisp.pl?gene="+geneinfo.genesym)

  root.find("#genepanel-toset").click () -> addtogeneset geneinfo.genesym

  ########### count histogram ##########

  ##Pull out data, log2-transform
  datasetlist = Object.keys(geneinfo["dataset"])
  values={}
  for k in datasetlist
    values[k] = geneinfo["dataset"][k]["exp"]
  values=assmap values, log2

  margin = {top: 10, right: 30, bottom: 30, left: 30}
  width = 650 - margin.left - margin.right
  height = 200 - margin.top - margin.bottom

  ##Find x-range
  maxx = d3.max (d3.values values), (arr)-> 
    d3.max arr, ((d)->d)
  minx = -2
  maxx = Math.max maxx,5
  x = d3.scale.linear()
    .domain([minx, maxx])   
    .range([0, width])
  #kde = kernelDensityEstimator(epanechnikovKernel(0.5), x.ticks(100))
  kde = kernelDensityEstimator(gaussianKernel(0.5), x.ticks(100))

  ##Create bins
  data={}
  for k in Object.keys(values)
    data[k]=kde(values[k])

  ##Compute y-range
  miny=0.0
  maxy=0.0
  for k in Object.keys(values)
    maxy=Math.max(maxy,d3.max data[k], (d)->d[1])
  y = d3.scale.linear()
    .domain([miny, maxy])
    .range([height, 0])

  xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom")
  yAxis = d3.svg.axis()
    .scale(y)
    .orient("left")

  svg = d3.select("#expressionhist")
  svg = svg
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")

  ##### draw all histograms
  thelines = {}
  for k in Object.keys(values)
    thelines[k] = oneline svg,x,y,data[k],thecol[k]
    if !(defaultshowds[k] == 1)
      thelines[k].style "opacity", 0

  svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis);

  svg.append("text")
    .attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom) + ")")
    .style("text-anchor", "middle")
    .text("Log normalized expression level");
  svg.append("text")
    .attr("transform", "rotate(-90)")
    .attr("y", 0 - margin.left)
    .attr("x",0 - (height / 2))
    .attr("dy", "1em")
    .style("text-anchor", "middle")
    .text("Density (# cells)");

  ##Add legend
  elegend = (root.find("#countlegend")) 
  for k in Object.keys(values).sort()
    etr = $ "<tr>"
    elegend.append etr
    etd = $ "<td>"
    etr.append etd

    echeck = $ "<input type=\"checkbox\">"
    if defaultshowds[k] == 1
      echeck.attr checked: "checked"
    etd.append echeck
    efont = $ "<font>"
    efont.attr "color", thecol[k]
    efont.append (" "+k)
    etd.append efont

    etd = $ "<td>"
    etr.append etd
    thedm = geneinfo.genedm[k]
    if thedm != null
      etd.append thedm
    else
      etd.append "N/A"

    echeck.click (()->
      thisk=k
      thischeck=echeck
      thisline=thelines[k]
      return () ->
        if thischeck.is(':checked')
          newop=1
        else
          newop=0
        thisline.style "opacity", newop)()


  ########### highest correlations ##########

  ## Fill list of correlation dataset
  ecorrdataset = form2.find("#corrdataset")
  for v in Object.keys(geneinfo["corr"])
    eoption = $ "<option>"
    eoption.append v
    eoption.attr "value", v
    ecorrdataset.append eoption    

  ##Fill table of correlations
  fillcorrtable = () ->
    ecorrdataset = form2.find("#corrdataset")
    values=geneinfo["corr"][ecorrdataset.val()]
    shownumcorr = +form2.find("#corrshownum").val()
    table = form2.find("#corrlist")
    table.empty()
    maxloop = 1
    if (typeof values) == "undefined"
      row = $ "<tr>"
      row.append "(not computed for this gene)"
      table.append row
    else
      values = values.sort (a,b) -> 
        Math.abs(b.corr)-Math.abs(a.corr)
      for rec in values
        row = $ "<tr>"
        row.append (onetd (makegenelink rec["geneid"], rec["geneid"],  true))
        row.append (onetd (makegenelink rec["geneid"], rec["genesym"], true))
        row.append (onetd rec["corr"])
        table.append row
        maxloop++
        if maxloop == shownumcorr
          break
  
  fillcorrtable()
  
  form2.find("#corrshownum").change () ->
    fillcorrtable()
  form2.find("#corrdataset").change () ->
    fillcorrtable()


oneline = (svg,x,y,data,color) ->
  line = d3.svg.line()
     .x((d) -> x(d[0]))
     .y((d) -> y(d[1]))
     #.x((d) -> x(d.x)) #for histogram data
     #.y((d) -> y(d.y))
  return svg.append("path")
    .datum(data)
    .attr("class","line")
    .attr("stroke-width","3px")
    .attr("fill","none")
    .attr("stroke",color)
    .attr("d", line)






###########################################################################
###########################################################################

assmap = (values, func) ->
  newvalues={}
  for k in Object.keys(values)
    newvalues[k]=values[k].map func
  return newvalues
  
onetd = (content) ->
  td = $ "<td class=\"corrtd\"/>"
  td.append content
  





kernelDensityEstimator = (kernel, x) ->
  (sample) ->
    x.map (x) ->
      [x, d3.mean sample, ((v) -> kernel (x - v))]

epanechnikovKernel = (scale) ->
  return (u) ->
    if Math.abs(u /= scale) <= 1
      return 0.75 * (1 - u * u) / scale
    else
      return 0

gaussianKernel = (sigma) ->
  thediv = 1.0/(sigma*Math.sqrt(2*Math.PI))
  thediv2 = 1.0/(2*sigma*sigma)
  return (u) ->
    return thediv*Math.exp(-u*u*thediv2)





###########################################################################
# Open up a gene set view
###########################################################################

add_genelist = (genelist) ->
  req = $.get genelist
  req.success (data) ->
    genesetlist = $ "#genesetlist"
    genesetlist.val (genesetlist.val()+"\n"+data).trim()
  req.fail (data) -> 
    alert("failed to query data, "+genelist)


onlyUnique = (value,index,self) ->
  (self.indexOf value) == index


getleftgeneset = () ->
  v = ($ "#genesetlist").val()
  v = (v.split "\n").join " "
  v = (v.split "\r").join " "
  v = (v.split "\t").join " "
  v = v.split " "
  v = v.filter ((d)->d!="")
  v = v.filter onlyUnique
  setquery =
    geneset: v



view_geneset = () ->
  usetab = current_divname
  closecurrentdiv()
  hide_startpage()
  setquery = getleftgeneset()

  if setquery.geneset.length<2 || setquery.geneset.length>50
    alert("Need to select more than one gene, and at the very most 50")
  else
    query_url = "geneset.php?q=#{JSON.stringify setquery}"   #todo, use POST. and different query!
    req = $.getJSON query_url
    req.success (data) ->
      view_geneset_disp(data,usetab)
    req.fail (data) -> 
      alert("failed to query data, "+query_url)
  


view_geneset_disp = (data, curtab) ->
  dispidcounter++
  thisid="genesetdisp"+dispidcounter

  root = $ "#maincontent"
  form2 = $ "#genesetpanel"
  form2 = form2.html()
  root.prepend form2
  form2 = root.children(":first")
  form2.id = thisid
  
  root.find("#genesetpanel-close").click () ->
    form2.remove()

  $("html, body").animate (scrollTop: form2.offset().top), "slow"


  updatetab = () ->
    form2.find("#tab_betweengenes").addClass "hideclass"
    form2.find("#tab_cellcorr").addClass "hideclass"
    form2.find("#btab_betweengenes").css  "background-color", tabcolorInactive
    form2.find("#btab_cellcorr").css "background-color", tabcolorInactive
    if(curtab=="cellcor")
      form2.find("#btab_cellcorr").css "background-color", tabcolorActive
      form2.find("#tab_cellcorr").removeClass "hideclass"
    else
      form2.find("#btab_betweengenes").css "background-color", tabcolorActive
      form2.find("#tab_betweengenes").removeClass "hideclass"
  updatetab()


  form2.find("#btab_cellcorr").click () ->
    curtab="cellcor"
    updatetab()
  form2.find("#btab_betweengenes").click () ->
    curtab="betweengenes"
    updatetab()



  getthisgeneset = () ->
    q = []
    for k in data.geneset
      q.push k.geneid
    return geneset: q

  egenegeneds = form2.find("#genegenedataset")
  for k in Object.keys(thecol)
    egenegeneds.append makeoption k,k


  ############# gene-gene correlation map
  egenesetcorr1 = form2.find("#genesetcorr1")
  updatecorr1 = () ->
    setquery = getthisgeneset() 
    setquery.graphw=2000
    setquery.dataset=form2.find("#genegenedataset").val() 
    query_url = "corrgenegene.php?q=#{JSON.stringify setquery}"
    #alert(query_url)
    imglink = $ "<a>"
    imglink.attr "href", query_url
    img = $ "<img>"

    setquery = getthisgeneset()
    setquery.dataset=form2.find("#genegenedataset").val() 
    query_url = "corrgenegene.php?q=#{JSON.stringify setquery}"
    img.attr src: query_url
    img.attr width:  700
    img.attr height: 700
    imglink.append img
    egenesetcorr1.empty()
    egenesetcorr1.append imglink

  updatecorr1()
  egenegeneds.change () ->
    updatecorr1()

  ############# gene-cell correlation map
  imglink = $ "<a>"
  img = $ "<img style=\"margin-right: -100; overflow: hidden;\">"
  img.attr width:  700
  img.attr height: 700
  imglink.append img
  root.find("#genesetcorr2").append imglink


  elegend = (root.find("#genesetcorr2legend")) 
  echecks = {}
  echecksid = {}
  for k in Object.keys(thecol)
    echecks[k] = $ "<input type=\"checkbox\">"
    if defaultshowds[k] == 1
      echecks[k].attr checked: "checked"
    dispidcounter++
    echecks[k].attr "id", "genedisp"+dispidcounter

    efont = $ "<font>"
    efont.attr "color", thecol[k]
    efont.append (" "+k)
    elegend.append echecks[k]
    elegend.append efont
    elegend.append "<br/>"

    echecksid[k] = echecks[k].attr("id")

    
  get_genecorr_queryurl = (()->
    thisimg = img
    thischecksid=echecksid
    setquery = getthisgeneset()
    return (graphw)->
      setquery.datasets = []
      for k in Object.keys(thischecksid)
        thischeck = $("#"+thischecksid[k])
        if thischeck.is(':checked')
          setquery.datasets.push(k)
      setquery.graphw=graphw
      query_url = "corrgenecell.php?q=#{JSON.stringify setquery}"
      return query_url
      )()

  update_genecorr = (()->
    thisimg = img
    thisimglink = imglink
    thischecksid=echecksid
    setquery = getthisgeneset()
    return ()->
      thisimglink.attr href: get_genecorr_queryurl(2000)
      thisimg.attr src: get_genecorr_queryurl(500)
      )()
  update_genecorr()

  for k in Object.keys(thecol)
    echeck = echecks[k]
    echeck.click (()->
      return () ->
        update_genecorr())()



  values=data["geneset"]
  egenesetlist = root.find("#genesetlist")
  egenesetlist.append "<b>"
  for rec in values
    egenesetlist.append (makegenelink rec.geneid, rec.genesym, true)
    egenesetlist.append " "





###########################################################################
# Initiate page
###########################################################################

$ ->
  $("#searchform").submit search_genes_form


  $("#search").keyup () ->
    query.gene = $(this).val()
    search_genes()

  search_genes()

  $('#view-diffexp').click () ->
    view_diffexp()
  $('#view-gonoise').click () ->
    view_godm()

  $('#view-geneset').click () ->
    view_geneset()

  $('#reset-geneset').click ()->
    genesetlist = $ "#genesetlist"
    genesetlist.val ""

  $('#genelist_pluripotency').click () ->
    add_genelist "genelist_pluripotency.txt"

  floatdivhandler '#openfloat-geneexp','#float-pickonegene',"geneexp"
  floatdivhandler '#openfloat-genecorrone','#float-pickonegene',"genecorrone"
  floatdivhandler '#openfloat-genecorrmany','#float-pickgeneset',GENESETTAB_BETWEENGENES
  floatdivhandler '#openfloat-cellcorr','#float-pickgeneset',GENESETTAB_CELLCOR
  floatdivhandler '#openfloat-download','#float-download',"download"

  $(document).mousemove (e) ->
   mouseX = e.pageX 
   mouseY = e.pageY

  $("#paperlink").attr href: publication_url

  reqsvg = $.get "startsite.svg"
  reqsvg.success (data) ->
    thesvg = $ "#svg2"
    thesvg.prepend (new XMLSerializer().serializeToString(data))

    floatdivhandler2 "#svg_geneexp","#float-pickonegene","geneexp"
    floatdivhandler2 "#svg_genecorrone","#float-pickonegene","genecorrone"
    floatdivhandler2 "#svg_genecorrmany","#float-pickgeneset",GENESETTAB_BETWEENGENES
    floatdivhandler2 "#svg_cellcorr","#float-pickgeneset",GENESETTAB_CELLCOR
    floatdivhandler2 "#svg_download","#float-download","download"

    ($ "#svg2").find("#svg_diffexp").click () ->
      view_diffexp()
    ($ "#svg2").find("#svg_publication").click () ->
      view_publication()
    ($ "#svg2").find("#svg_gonoise").click () ->
      view_godm()


view_publication = () ->
  window.location.href = publication_url

hide_startpage = () ->
  $("#startpage").addClass "hideclass"
  $("#topmenu").removeClass "hideclass"


getMethods = (obj) ->
  result = []
  for id in obj
    #if (typeof(obj[id]) == "function")
    result.push(id + ": " + obj[id].toString())
  return result



##### Set a button to open a hidden div
floatdivhandler = (buttonid,floatname,divname) ->
  thebutton = $(buttonid)
  thebutton.click () ->
    thediv = $(floatname)

    if(current_divname!=divname)
      closecurrentdiv()

    if(thediv.hasClass "hideclass")
      thediv.removeClass "hideclass"
      thediv.css "top",(thebutton.position().top+thebutton.outerHeight())
      thediv.css "left",(thebutton.position().left)
      current_div=thediv
      current_divname=divname
    else
      thediv.addClass "hideclass"
      current_div=null
      current_divname=divname


##### Set a button to open a hidden div   (called from link in div)
floatdivhandler2 = (linkid, floatname,divname) ->
  thesvg = $ "#svg2"
  thesvg.find(linkid).click () ->
    thediv = $(floatname)
    if(current_divname!=divname)
      closecurrentdiv()

    if(thediv.hasClass "hideclass")
      thediv.removeClass "hideclass"
      thediv.css "top",mouseY
      thediv.css "left",mouseX
      current_div=thediv
      current_divname=divname
    else
      thediv.addClass "hideclass"
      current_div=null
      current_divname=divname




#### Close current div
closecurrentdiv = () ->
  if(current_div!=null)
    current_div.addClass "hideclass"
    current_div=null
    current_divname=null




















###########################################################################
# Show differential expression
###########################################################################


view_diffexp = () ->
  closecurrentdiv()
  hide_startpage()
  theurl = "availablediffexp.php"
  req = $.getJSON theurl
  req.success (data) ->
    view_diffexp2 data
  req.fail (data) -> 
    alert("failed to query data, "+theurl)


view_diffexp2 = (dslist) ->
  dispidcounter++
  thisid="genedisp"+dispidcounter

  root = $ "#maincontent"
  form2 = $ "#diffexppanel"
  form2 = form2.html()
  root.prepend form2
  form2 = root.children(":first")
  form2.attr id: thisid
  
  $("html, body").animate (scrollTop: form2.offset().top), "slow"

  root.find("#diffexppanel-close").click () ->
    form2.remove()

  eds1 = form2.find("#ds1")
  eds2 = form2.find("#ds2")
  etable = form2.find("#diffexptable")
  egraph = form2.find("#diffexpgraph")

  getdsquery = () ->
    query =
      dataset1: eds1.val()
      dataset2: eds2.val()
    return "diffexp.php?q=#{JSON.stringify query}"   

  elinkcsv = form2.find("#linkcsv")
  updatecsvlink = () ->
    elinkcsv.attr href: getdsquery()

  ## Choice of dataset 1 and 2
  eds1.empty()
  for rec in (makeunique (dslist.map (x) -> x.dataset1)).sort()
    eds1.append makeoption rec,rec
  updateds2 = () ->
    eds2.empty()
    for rec in (makeunique ((dslist.filter (x) -> eds1.val()==x.dataset1).map (x) -> x.dataset2)).sort()
      eds2.append makeoption rec,rec
  updateds2()
  form2.find("#ds1").change () ->
    updateds2()

  ## Update the contents of the table
  filltable = (() ->
    thisform2 = form2
    return () ->
      shownumgenes = +thisform2.find("#diffexpshownum").val()
      d3.csv getdsquery(), (data) ->
        etable.empty()
        for i in [0..(Math.min shownumgenes,(data.length-1))]
          (() ->
            rec = data[i]
            rec1 = rec
      
            ea = $ "<a href=\"#\">(+)</a>"
            ea.click () ->
              addtogeneset rec1.genesym
              return false
            etd = $ "<td class=\"corrtd\"/>"
            etd.append ea
            etd.append " "
            etd.append (makegenelink rec.geneid, rec.genesym, true)
  
            etr = $ "<tr/>"
            etr.append etd
            etr.append onetd (log2 rec.mean1).toPrecision 2
            etr.append onetd (log2 rec.mean2).toPrecision 2
            etr.append onetd (log2 rec.mean2/rec.mean1).toPrecision 2
      
            etr.append onetd (+rec.pvalue).toPrecision 2
            etr.append onetd (+rec.padj).toPrecision 2
            etable.append etr
            )()


        ##Update the graph
        egraph.empty()

        margin = {top: 10, right: 30, bottom: 30, left: 30}
        width = 800 - margin.left - margin.right
        height = 400 - margin.top - margin.bottom
        svg = d3.select("#diffexpgraph")
        svg = svg
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
          .append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")")

        datasub = data[0..(Math.min shownumgenes,(data.length-1))]
        datasub = datasub.map (d) ->
          geneid: d.geneid
          genesym: d.genesym
          mean1: Math.max -3,(log2 d.mean1)
          mean2: Math.max -3,(log2 d.mean2)

        minx = -1 + (d3.min datasub, ((d)->Math.max -3,d.mean1))
        miny = -1 + (d3.min datasub, ((d)->Math.max -3,d.mean2))
        maxx =  1 + (d3.max datasub, ((d)->d.mean1))
        maxy =  1 + (d3.max datasub, ((d)->d.mean2))
        x = d3.scale.linear()
          .domain([minx, maxx])   
          .range([0, width])
        y = d3.scale.linear()
          .domain([miny, maxy])
          .range([height, 0])
        xAxis = d3.svg.axis()
          .scale(x)
          .orient("bottom")
        yAxis = d3.svg.axis()
          .scale(y)
          .orient("left")
        svg.append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(0," + (height) + ")")
          .call(xAxis);
        svg.append("g")
          .attr("transform", "translate(5,0)")
          .attr("class", "y axis")
          .call(yAxis);
        svg.append("text")
          .attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom) + ")")
          .attr("font-size","10px")
          .style("text-anchor", "middle")
          .text("Log mean (from)");
        svg.append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", 0 - margin.left)
          .attr("x",0 - (height / 2))
          .attr("dy", "1em")
          .attr("font-size","10px")
          .style("text-anchor", "middle")
          .text("Log mean (to)");

        for i in [0..(datasub.length-1)]
          svg.append("circle")
            .attr("cx", x datasub[i].mean1)
            .attr("cy", y datasub[i].mean2)
            .attr("r", 3)
            .style("fill", "red")
          ea = svg.append("a")
            .on("click",makeviewgenefunc(datasub[i].geneid),true)
            .append("text")
            .attr("x", x datasub[i].mean1)
            .attr("y", (y datasub[i].mean2)-5)
            .style("font-weight","bold")
            .style("text-anchor","middle")
            .text(datasub[i].genesym)
      )()

      

  filltable()
  updatecsvlink()

  form2.find("#diffexpshownum").change () ->
    filltable()
  form2.find("#ds1").change () ->
    filltable()
    updatecsvlink()
  form2.find("#ds2").change () ->
    filltable()
    updatecsvlink()



makeoption = (val,title) ->
  eopt = $ "<option/>"
  eopt.attr value: val
  eopt.append title
  return eopt


makeunique = (a) ->
  temp = {}
  for x in a
    temp[x] = true
  r = []
  for k in Object.keys(temp)
    r.push k
  return r


###########################################################################
# Show GO DM
###########################################################################


view_godm = () ->
  closecurrentdiv()
  hide_startpage()
  theurl = "availablegodm.php"
  req = $.getJSON theurl
  req.success (data) ->
    view_godm2 data
  req.fail (data) -> 
    alert("failed to query data, "+theurl)


view_godm2 = (dslist) ->
  dispidcounter++
  thisid="genedisp"+dispidcounter

  root = $ "#maincontent"
  form2 = $ "#godmpanel"
  form2 = form2.html()
  root.prepend form2
  form2 = root.children(":first")
  form2.attr id: thisid
  
  $("html, body").animate (scrollTop: form2.offset().top), "slow"

  form2.find("#bclose").click () ->
    form2.remove()

  eds1 = form2.find("#ds1")
  eds2 = form2.find("#ds2")
  etable = form2.find("#godmtable")

  getdsquery = () ->
    query =
      dataset1: eds1.val()
      dataset2: eds2.val()
    return "godm.php?q=#{JSON.stringify query}"   

  elinkcsv = form2.find("#linkcsv")
  updatecsvlink = () ->
    elinkcsv.attr href: getdsquery()

  ## Choice of dataset 1 and 2
  eds1.empty()
  for rec in (makeunique (dslist.map (x) -> x.dataset1)).sort()
    eds1.append makeoption rec,rec
  updateds2 = () ->
    eds2.empty()
    for rec in (makeunique ((dslist.filter (x) -> eds1.val()==x.dataset1).map (x) -> x.dataset2)).sort()
      eds2.append makeoption rec,rec
  updateds2()
  form2.find("#ds1").change () ->
    updateds2()

  ## Update the contents of the table
  filltable = (() ->
    thisform2 = form2
    return () ->
      shownum = +thisform2.find("#godmshownum").val()
      d3.csv getdsquery(), (data) ->
        etable.empty()
        for i in [0..(Math.min shownum,(data.length-1))]
          (() ->
            rec = data[i]
            etr = $ "<tr/>"
            etr.append onetd rec.goid
            etr.append onetd rec.goname
            etr.append onetd rec.pvalue
            etr.append onetd (+rec.tscore).toPrecision 2
            etable.append etr
            )()
      )()

      

  filltable()
  updatecsvlink()

  form2.find("#godmshownum").change () ->
    filltable()
  form2.find("#ds1").change () ->
    filltable()
    updatecsvlink()
  form2.find("#ds2").change () ->
    filltable()
    updatecsvlink()
