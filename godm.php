<?php
include "common.php";

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$ds1 = filteralpha(getdef($query, 'dataset1', ''));
$ds2 = filteralpha(getdef($query, 'dataset2', ''));


if($ds1=="" && $ds2==""){
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM godm NATURAL JOIN goinfo");
  $rs=pg_execute($db, 'getgene', array());
} else {
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM godm NATURAL JOIN goinfo WHERE dataset1=$1 AND dataset2=$2 ORDER BY pvalue");
  $rs=pg_execute($db, 'getgene', array($ds1,$ds2));
}
$results=array();
echo '"dataset1","dataset2","goid","goname","pvalue","tscore"'."\n";
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  echo '"' . $line["dataset1"] . '","' . $line["dataset2"] . '","' . 
    $line["goid"] . '","' . $line["goname"] . '","' . $line["pvalue"] . '","' . $line["tscore"] . "\"\n";
}
pg_free_result($rs);



### Return result
header('Content-Type: application/json; ');
header('filename="godm.json"; ');
echo json_encode($results);


include "end.php";
?>
