<?php
include "common.php";

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$ds1 = filteralpha(getdef($query, 'dataset1', ''));
$ds2 = filteralpha(getdef($query, 'dataset2', ''));

echo "dataset1,dataset2,geneid,genesym,mean1,mean2,pvalue,padj\n";

if($ds1=="" && $ds2==""){
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM diffexp NATURAL JOIN (SELECT geneid,genesym FROM geneinfo) as bar");
  $rs=pg_execute($db, 'getgene', array());
} else {
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM diffexp NATURAL JOIN (SELECT geneid,genesym FROM geneinfo) as bar WHERE dataset1=$1 AND dataset2=$2 ORDER BY padj");
  $rs=pg_execute($db, 'getgene', array($ds1,$ds2));
}
$results=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  echo '"' . $line["dataset1"] . '","' . $line["dataset2"] . '","' . 
        $line["geneid"] . '","' . $line["genesym"] . '","' . $line["mean1"] . '","' . $line["mean2"] . '","' . 
        $line["pvalue"] . '","' . $line["padj"] . "\"\n";
}
pg_free_result($rs);



### Return result
header('Content-Type: application/json; ');
echo json_encode($results);


include "end.php";
?>
