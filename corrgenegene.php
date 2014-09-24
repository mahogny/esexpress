<?php

  include "common.php";

  $query = json_decode(getdef($_GET,'q','{}'),TRUE);
  $geneid = getdef($query, 'geneset', array("ENSMUSG00000000126","ENSMUSG00000000028"));
  $graphw = getdef($query, 'graphw', 500);
  $dataset = getdef($query, 'dataset', 'es_lif')

#freaking dangerous command!######## TODO clean

  $cmd = "echo 'genes <- c();";
  foreach($geneid as $s){
    $cmd=$cmd . "genes <- c(genes,\"".$s."\");";
  }
  $cmd = $cmd . "graphw<-" . $graphw . ";";
  $cmd = $cmd . "dataset<-" . $dataset . ";";
  $cmd = $cmd . "source(\"corrgenegene.R\")' | /usr/bin/R --vanilla --slave";

  #echo $cmd;

  $handle = popen($cmd, "r");
  $ret = "";
  do{
    $data = fread($handle, 8192);
    if(strlen($data) == 0){
      break;
    }
    $ret .= $data;
  }
  while(true);
  pclose($handle);

  header("Content-type:image/png");
  header("Cache-Control: no-cache, must-revalidate"); 
  echo $ret;

?>
