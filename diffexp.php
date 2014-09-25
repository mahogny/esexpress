<?php
include "common.php";
header('Content-Type: application/json; ');

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$ds1 = filteralpha(getdef($query, 'dataset1', ''));
$ds2 = filteralpha(getdef($query, 'dataset2', ''));

echo "geneid,genesym,mean1,mean2,pvalue,padj\n";

$ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM diffexp NATURAL JOIN (SELECT geneid,genesym FROM geneinfo) as bar WHERE dataset1=$1 AND dataset2=$2 ORDER BY padj");
$rs=pg_execute($db, 'getgene', array($ds1,$ds2));
$results=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  echo $line["geneid"] . ',' . $line["genesym"] 
        . ',' . $line["mean1"] . ',' . $line["mean2"] . ',' . $line["pvalue"] . ',' . $line["padj"] . "\n";
}
pg_free_result($rs);



### Return result
echo json_encode($results);


include "end.php";
?>
