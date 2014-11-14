<?php
include "common.php";

$query = json_decode(getdef($_GET,'q','{}'),TRUE);
$ds1 = filteralpha(getdef($query, 'dataset', ''));

if($ds1==""){
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM genedm");
  $rs=pg_execute($db, 'getgene', array());
} else {
  $ps=pg_prepare($db, 'getgene',"SELECT DISTINCT * FROM genedm WHERE dataset=$1");
  $rs=pg_execute($db, 'getgene', array($ds1));
}
echo '"dataset","geneid","genedm"'."\n";
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC)){
  echo '"' . $line["dataset"] . '","' . $line["geneid"] . '","' .  $line["genedm"] . "\"\n";
}
pg_free_result($rs);



### Return result
#header('Content-Type: application/json; ');
#header('filename="godm.json"; ');


include "end.php";
?>
