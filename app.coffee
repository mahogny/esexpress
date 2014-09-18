#datasets = ["2i", "a2i", "lif", "s_lblast", "s_mblast", "s_eblast"]

query = 
  gene: ""
dispidcounter=1

thecol = {}
thecol["es_lif"] = "#8B0000"
thecol["es_a2i"] = "#EEAD0E"
thecol["es_2i"]  = "#00008B"

thecol["sandberg_earlyblast"] = "#FF0000"
thecol["sandberg_midblast"]   = "#00FF00"
thecol["sandberg_lateblast"]  = "#0000FF"


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
    link = makegenelink rec.geneid, rec.genesym
    table.append link
    table.append ($ "<br/>")



makegenelink = (geneid, genesym) ->
    link = $ "<a href=\"#\">"
    link.append genesym
    link.click (() -> 
      return () ->
        view_gene(geneid)
        #$("html, body").animate (scrollTop: 0), "slow"
        return false)()








###########################################################################
# Open up a single gene view
###########################################################################

view_gene = (gene) ->
  query.gene = gene
  query_url = "geneinfo.php?q=#{JSON.stringify query}"   
  req = $.getJSON query_url
  req.success (data) ->
    view_gene_disp(data)
  req.fail (data) -> 
    alert("failed to query data")



view_gene_disp = (geneinfo) ->
  dispidcounter++
  thisid="genedisp"+dispidcounter

  root = $ "#maincontent"
  form2 = $ "#genepanel"
  form2 = form2.html()
  form2 = form2.replace("GENEID", geneinfo.geneid)
  form2 = form2.replace("GENESYM", geneinfo.genesym)
  root.prepend form2
  form2 = root.children(":first")
  #  form2.id = thisid
  form2.attr id: thisid
  
  $("html, body").animate (scrollTop: form2.offset().top), "slow"

  root.find("#genepanel-close").click () ->
    form2.remove()

  form2.find("#linkpubmed").attr href: ("http://www.ncbi.nlm.nih.gov/pubmed/?term="+geneinfo.genesym)
  form2.find("#linkensembl").attr href: ("http://www.ensembl.org/Mus_musculus/Gene/Summary?g="+geneinfo.geneid)
  form2.find("#linkgenecards").attr href: ("http://www.genecards.org/cgi-bin/carddisp.pl?gene="+geneinfo.genesym)

  root.find("#genepanel-toset").click () ->
    genesetlist = $ "#genesetlist"
    genesetlist.val (genesetlist.val()+"\n"+geneinfo.genesym).trim()

  ########### count histogram ##########

  ##Pull out data, log2-transform
  datasetlist = Object.keys(geneinfo["dataset"])
  values={}
  for k in datasetlist
    values[k] = geneinfo["dataset"][k]["exp"]
  values=assmap values, Math.log2

  margin = {top: 10, right: 30, bottom: 30, left: 30}
  width = 400 - margin.left - margin.right
  height = 200 - margin.top - margin.bottom

  ##Find x-range
  maxx = d3.max (d3.values values), (arr)-> 
    d3.max arr, ((d)->d)
  minx = -3
  #maxx = 11
  x = d3.scale.linear()
    .domain([minx, maxx])   
    .range([0, width])

  ##Create bins
  data={}
  for k in Object.keys(values)
    data[k]=d3.layout.histogram()
              .bins(x.ticks(30))(values[k])

  ##Compute y-range
  miny=0
  maxy=0
  for k in Object.keys(values)
    maxy=Math.max(maxy,d3.max data[k], (d)->d.y)
  y = d3.scale.linear()
    .domain([miny, maxy])
    .range([height, 0])

  xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom")

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

  svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis);

  ##Add legend
  elegend = (root.find("#countlegend")) 
  for k in Object.keys(values)
    echeck = $ "<input type=\"checkbox\" checked=\"checked\">"
    elegend.append echeck
    efont = $ "<font>"
    efont.attr "color", thecol[k]
    efont.append (" "+k)
    elegend.append efont
    elegend.append "<br/>"

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
    #alert(JSON.stringify values)
    if (typeof values) == "undefined"
      row = $ "<tr>"
      row.append "(not computed)"
      table.append row
    else
      for rec in values
        row = $ "<tr>"
        row.append (onetd (makegenelink rec["geneid"], rec["genesym"]))
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
     .x((d) -> x(d.x))
     .y((d) -> y(d.y))
  return svg.append("path")
    .datum(data)
    .attr("class","line")
    .attr("stroke-width","1px")
    .attr("fill","none")
    .attr("stroke",color)
    .attr("d", line)



assmap = (values, func) ->
  newvalues={}
  for k in Object.keys(values)
    newvalues[k]=values[k].map func
  return newvalues
  


onetd = (content) ->
  td = $ "<td class=\"corrtd\"/>"
  td.append content
  




###########################################################################
# Open up a gene set view
###########################################################################

add_genelist = (genelist) ->
    req = $.get genelist
    req.success (data) ->
      genesetlist = $ "#genesetlist"
      genesetlist.val (genesetlist.val()+"\n"+data).trim()
    req.fail (data) -> 
      alert("failed to query data")


getgeneset = () ->
  v = ($ "#genesetlist").val()
  v = (v.split "\n").join " "
  v = (v.split "\r").join " "
  v = (v.split "\t").join " "
  v = v.split " "
  v = v.filter ((d)->d!="")
  setquery =
    geneset: v

view_geneset = () ->
  setquery = getgeneset()
  if setquery.geneset.length<2
    alert("Need to select more than one gene")
  else
    query_url = "geneset.php?q=#{JSON.stringify setquery}"   #todo, use POST. and different query!
    req = $.getJSON query_url
    req.success (data) ->
      view_geneset_disp(data)
    req.fail (data) -> 
      alert("failed to query data, "+query_url)
  


view_geneset_disp = (data) ->

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

  setquery = getgeneset()
  query_url = "corrgenegene.php?q=#{JSON.stringify setquery}"   #todo, use POST. and different query!
  img = $ "<img>"
  img.attr src: query_url
  img.attr width:  500
  img.attr height: 500
  root.find("#genesetcorr1").append img

  ############# gene-cell correlation map
  imglink = $ "<a>"
  img = $ "<img style=\"margin-right: -100; overflow: hidden;\">"
  img.attr width:  500
  img.attr height: 500
  imglink.append img
  root.find("#genesetcorr2").append imglink


  elegend = (root.find("#genesetcorr2")) 
  echecks = {}
  echecksid = {}
  for k in Object.keys(thecol)
    echecks[k] = $ "<input type=\"checkbox\" checked=\"checked\">"
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
    setquery = getgeneset()
    return (graphw)->
      setquery.datasets = []
      for k in Object.keys(thischecksid)
        thischeck = $("#"+thischecksid[k])
        if thischeck.is(':checked')
          setquery.datasets.push(k)
      setquery.graphw=graphw
      query_url = "corrgenecell.php?q=#{JSON.stringify setquery}"   #todo, use POST. and different query!
      return query_url
      )()

  update_genecorr = (()->
    thisimg = img
    thisimglink = imglink
    thischecksid=echecksid
    setquery = getgeneset()
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
    egenesetlist.append (makegenelink rec.geneid, rec.genesym)
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

  $('#view-geneset').click view_geneset
  $('#reset-geneset').click ()->
    genesetlist = $ "#genesetlist"
    genesetlist.val ""

  $('#genelist_pluripotency').click () ->
    add_genelist "genelist_pluripotency.txt"






