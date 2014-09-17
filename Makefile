all:
	coffee -o static/ -c .
clean:
	rm static/app.js

#upload:
#	scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" -r * ebi-001.ebi.ac.uk:/nfs/public/rw/research/teichmann/thexpress/www/


