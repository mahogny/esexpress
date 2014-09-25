<?php

## http://localhost/geneinfo.php?q={%22gene%22:%22ENSMUSG00000000049%22}

include "common.php";

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$geneid = filteralpha(getdef($query, 'gene', ''));


### Query database for expression counts, for all datasets
$time_getgene = microtime();
$ps=pg_prepare($db, 'getgene','SELECT * FROM geneexp WHERE fromgene=$1');
$rs=pg_execute($db, 'getgene', array($geneid));
$resultsexp=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $resultsexp[$line['dataset']]=array('fromcell' => splitcomma($line['fromcell']), 'exp' => splitcomma($line['exp']));
        }
pg_free_result($rs);
$time_getgene = microtime()-$time_getgene;

$results=array();
$results['dataset']=$resultsexp;
$results['geneid']=$geneid;


### Query database for gene info --- this particular gene
$time_getgenei = microtime();
$ps=pg_prepare($db, 'getgenei','SELECT * FROM geneinfo WHERE geneid=$1 LIMIT 1');
$rs=pg_execute($db, 'getgenei', array($geneid));
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $results['genesym']=$line['genesym'];
        }
pg_free_result($rs);
$time_getgenei = microtime()-$time_getgenei;


### Query database for a map ensemblid -> genesym
$time_getgeneall = microtime();
$ps=pg_prepare($db, 'getgeneall','SELECT geneid,genesym FROM geneinfo');
$rs=pg_execute($db, 'getgeneall', array());
$listofgenes=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
  $listofgenes[$line["geneid"]]=$line["genesym"];
pg_free_result($rs);
$time_getgenei = microtime()-$time_getgenei;


### Get all gene correlations
$time_getcorr = microtime();
pg_prepare($db, 'select1','SELECT * from genecorr WHERE fromgene=$1');
$rs=pg_execute($db, 'select1', array($geneid));
$outcorr=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
  {
  $corr=splitcomma($line['corr']);
  $togene=splitcomma($line['togene']);
  for($i=0;$i<sizeof($togene);$i++)
    {
    $togeneid = $togene[$i];
    if(array_key_exists($togeneid,$listofgenes))
      $togenesym = $listofgenes[$togeneid];
    else
      $togenesym = "missing in db ".$togeneid;
    $thiscorr = $corr[$i];
    if($togeneid != $geneid)
      $outcorr[$line['dataset']][]=array('geneid' => $togeneid, 'genesym' => $togenesym, 'corr' => $thiscorr);
    }
  }
pg_free_result($rs);
$results['corr']=$outcorr;
$time_getcorr = microtime() - $time_getcorr;


$results['time_getgene']=$time_getgene;
$results['time_getgenei']=$time_getgenei;
$results['time_getcorr']=$time_getcorr;




### Return result
header('filename="data.json"; ');
header('Content-Type: application/json; ');
echo json_encode($results);


include "end.php";


?>
