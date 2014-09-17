<?php

include "common.php";

### Query database
#$q = getdef($_GET,'q','');
#$gene = $q['gene'];
$gene = getdef($_GET, 'gene', '');

#echo $gene;

if($gene!=''){
#	$ps=pg_prepare($db, 'getgene','SELECT * FROM geneinfo WHERE to_tsvector(genesym || \' \' || geneid) @@ to_tsquery($1) LIMIT 20');
  $ps=pg_prepare($db, 'getgene','SELECT *, similarity(genesym,$1) as sim FROM geneinfo WHERE similarity(genesym,$1)>0.1 OR $1=geneid LIMIT 20');
	$rs=pg_execute($db, 'getgene', array($gene));
} else {
	$ps=pg_prepare($db, 'getgene','SELECT *, -1 as sim FROM geneinfo LIMIT 10');
	$rs=pg_execute($db, 'getgene', array());
}

### Fetch result, make JSON
$results=array();
while($line=pg_fetch_array($rs,null,PGSQL_ASSOC))
	{
	$results[]=array('geneid' => $line['geneid'], 'genesym' => $line['genesym'], 'sim' => $line['sim']);
	}
pg_free_result($rs);

### Return result
header('filename="data.json"; ');
header('Content-Type: application/json; ');
echo json_encode($results);  



include "end.php";
?>
