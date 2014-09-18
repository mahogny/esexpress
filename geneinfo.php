<?php

## http://localhost/geneinfo.php?q={%22gene%22:%22ENSMUSG00000000049%22}

include "common.php";


#then here, also query for the most correlated genes and slap the two together
#optionally, add a flag to geneinfo to include this or not





### Query database for expression counts, for all datasets
$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$geneid = getdef($query, 'gene', '');
$ps=pg_prepare($db, 'getgene','SELECT * FROM geneexp WHERE fromgene=$1');
$rs=pg_execute($db, 'getgene', array($geneid));
$resultsexp=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $resultsexp[$line['dataset']]=array('fromcell' => splitcomma($line['fromcell']), 'exp' => splitcomma($line['exp']));
        }
pg_free_result($rs);

$results=array();
$results['dataset']=$resultsexp;
$results['geneid']=$geneid;



### Query database for gene info
$ps=pg_prepare($db, 'getgenei','SELECT * FROM geneinfo WHERE geneid=$1 LIMIT 1');
$rs=pg_execute($db, 'getgenei', array($geneid));
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $results['genesym']=$line['genesym'];
        }
pg_free_result($rs);




### Query for most correlated genes
#$ps=pg_prepare($db, 'getgenecorr','SELECT * FROM genecorr WHERE fromgene=$1');
#$rs=pg_execute($db, 'getgenecorr', array($geneid));
#$resultscorr=array();
#while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
#        {
#        $togene=splitcomma($line['togene']);
#        $corr=splitcomma($line['corr']);
#        $arr=array();
#        for($i=0;$i<sizeof($togene);$i++)
#         {
#         $arr[]=array($togene[$i],$corr[$i]);
##$arr[$togene[$i]]=$corr[$i];
#         }
#        $resultscorr[$line['dataset']]=$arr;
# #       $resultscorr[$line['dataset']]=$arr;
#        }
#pg_free_result($rs);
#$results['corr']=$resultscorr;


#expand the array into a temporary table. is there a one-liner here?
$ps=pg_prepare($db, 'create1','create temporary table comparecorr(dataset TEXT, togene TEXT, corr REAL)');
pg_execute($db, 'create1', array());
pg_prepare($db, 'insert1','INSERT INTO comparecorr VALUES($1,$2,$3)');
pg_prepare($db, 'getgenecorr','SELECT * FROM genecorr WHERE fromgene=$1');
$rs=pg_execute($db, 'getgenecorr', array($geneid));
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $togene=splitcomma($line['togene']);
        $corr=splitcomma($line['corr']);
        for($i=0;$i<sizeof($togene);$i++)
          pg_execute($db, 'insert1', array($line['dataset'], $togene[$i], $corr[$i]));
        }
pg_free_result($rs);
pg_prepare($db, 'select1','SELECT * from comparecorr NATURAL JOIN (SELECT geneid AS togene, genesym FROM geneinfo) as foo ORDER BY abs(corr) DESC');
$rs=pg_execute($db, 'select1', array());
$outcorr=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
        {
        $togene = $line['togene'];
        if($togene != $geneid)
          $outcorr[$line['dataset']][]=array('geneid' => $togene, 'genesym' => $line['genesym'], 'corr' => $line['corr']);
        }
pg_free_result($rs);
$results['corr']=$outcorr;


#print($outcorr);



### Return result
header('filename="data.json"; ');
header('Content-Type: application/json; ');
echo json_encode($results);


#       echo prettyPrint($results);  #, indent=2, separators=(',', ': ')




include "end.php";


?>
