from rgt.Util import GenomeData
from rgt.helper import get_chrom_sizes_as_genomicregionset
g = GenomeData('hg38') 
regionset = get_chrom_sizes_as_genomicregionset(g.get_chromosome_sizes())

from rgt.ODIN.get_extension_size import get_extension_size
bamfile = 'PU1_cDC.chr19.bam'
ext, _ = get_extension_size(bamfile, start=0, end=300, stepsize=5)
