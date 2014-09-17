<?php

  include "common.php";

  $query = json_decode(getdef($_GET,'q','{}'),TRUE);
  $geneid = getdef($query, 'geneset', array("ENSMUSG00000000126","ENSMUSG00000000028"));


#freaking dangerous command!######## TODO clean

  $cmd = "echo 'genes <- c();";
  foreach($geneid as $s){
    $cmd=$cmd . "genes <- c(genes,\"".$s."\");";
  }
  
  $cmd = $cmd. "source(\"geneset.R\")' | /usr/bin/R --vanilla --slave";

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

  header('Content-Type: application/json; ');
  header("Cache-Control: no-cache, must-revalidate"); 
  echo $ret;

?>
