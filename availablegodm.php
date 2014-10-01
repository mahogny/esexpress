<?php
include "common.php";

$ps=pg_prepare($db, 'getgene',"SELECT DISTINCT dataset1,dataset2 FROM godm");
$rs=pg_execute($db, 'getgene', array());
$results=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  $results[]=$line;
}
pg_free_result($rs);



### Return result
#header('filename="data.json"; ');
header('Content-Type: application/json; ');
echo json_encode($results);


include "end.php";
?>
  