<?php

$myfile = fopen("auth.txt", "r") or die("Unable to open file!");
$dbname=fgets($myfile);
$dbuser=fgets($myfile);
$dbpass=fgets($myfile);
fclose($myfile);

$db = pg_connect("host=localhost dbname=$dbname user=$dbuser password=$dbpass") or die('could not connect to database');

###############################################################################
# Get a field or a default value
function getdef($col, $var,$def)
        {
        if(array_key_exists($var, $col))
                {
                $value = $col[$var];
                return $value;
                }
        else
                return $def;
        }



###############################################################################
# Read specified genes from a binary table file
function readbtable($tfile, $genelist, $numgenes)
	{
	$f = fopen($tfile,"rb");
	$out = array();
	for($i=0;$i<length($genelist);$i++)
		{
		fseek($f,$genelist[$i]*$numgenes-1); #from start???
		$out[]=fread($f,$numgenes); #double somehow
		}
	fclose($f);
	return $out;
	}


###############
# take comma-list to array, double
function splitcomma($x){
  $x = substr($x,1,strlen($x)-2);
	$x = explode(',',$x);
	return $x;
}


?>
