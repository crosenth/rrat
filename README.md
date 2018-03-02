This package contains one script that calculates median 16s copy numbers for
every tax_id available in NCBI (or --nodes equivalent csv).  Initial copy number
data is downloaded from the University of Michigain rrnDB project

https://rrndb.umms.med.umich.edu/static/download/

or equivalent passed as an argument. The NCBI taxonomy data is downloaded
from

ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

or passed as arguments --nodes and --merged in csv format.
The `merged` data translates old tax_ids that have been changed in 
latest nodes data.  The result is a csv file with columns tax_id,median.
