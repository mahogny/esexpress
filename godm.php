<?php
include "common.php";
header('Content-Type: application/json; ');

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$ds1 = filteralpha(getdef($query, 'dataset1', ''));
$ds2 = filteralpha(getdef($query, 'dataset2', ''));

echo '"goid","goname","pvalue","tscore"'."\n";

$ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM godm NATURAL JOIN goinfo WHERE dataset1=$1 AND dataset2=$2 ORDER BY pvalue");
$rs=pg_execute($db, 'getgene', array($ds1,$ds2));
$results=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  echo '"' . $line["goid"] . '","' . $line["goname"] . '","' . $line["pvalue"] . '","' . $line["tscore"] . "\"\n";
}
pg_free_result($rs);



### Return result
echo json_encode($results);


include "end.php";
?>
