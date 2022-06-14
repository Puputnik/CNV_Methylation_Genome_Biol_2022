# CNV_Methylation_Genome_Biol_2022

Requires R packages: "liftOver", "vroom", "gwascat", "tidyr", "ggplot2", "ggbeeswarm", "ggpubr", "RColorBrewer"

## PMD and Genome average methylation calculation in 10MB bins and in Copy Number Altered regions 

To download segmentation results, remora and deepsignal methylation raw files refere to Katsman et al publication.
The Genome_Biol_Data_analyzed folder contains the analysed data from Katsman et al.

Put hg19 segmentation results from NanoGladiator in ~/CNV_meth/data/CNV/ 
You can use any other software for CNV analysis as soon as you used hg19 as reference genome and the segmentation results files are tab separated files formatted in this way:

1) contig
2) start
3) end
4) number of bins for CNV calculation (not mandatory)
5) type of CNV (Amplification, Deletion, not mandatory)
6) predicted Copy Number (not mandatory)
7) Allelic Fraction (not mandatory)
8) Allelic Fraction error (not mandatory)
9) Log2ratio Segment Mean (MANDATORY)

```
1       890000  120450000       1142    Duplication     3       0.05    0.02    0.0191713516162252
1       145430000       204350000       546     Duplication     3       0.96    0.02    0.520735925624906
1       204450000       249110000       417     Duplication     3       0.4     0.03    0.238323946524685
2       60000   88955000        861     Deletion        1       0.24    0.02    -0.201281111783523
2       95630000        242900000       1394    Deletion        1       0.24    0.01    -0.211799130979085
...
Y       59017656        59217656        3       Deletion        1       0.52    0.03    -0.394951571943842
```
create a 2 column dict_chr.tsv file with:
1) Contig in the reference used for methylation analysis
2) Contig in chrN format

You can find the dict_chr.tsv files used for Katsman et al. in ~/CNV_meth/Utility/
For Deepsignal methylation analysis we used GCF_000001405.39_GRCh38.p13 as reference genome.
~/CNV_meth/Utility/dict_chr_DEEPSIGNAL.tsv
```
NC_000001.11	chr1
NC_000002.12	chr2
NC_000003.12	chr3
...
NC_000024.10	chrY
```
For Remora contigs were already in chrN format:
~/CNV_meth/Utility/dict_chr_REMORA.tsv
```
chr1	chr1
chr2	chr2
chr3	chr3
...
chrY	chrY
```

### Remora methylation data

Put hg38 remora methylation .bed files in ~/CNV_meth/data/REMORA/
```
chr1	10470	10471	5mC	1000	+	10470	10471	0,0,0	1	0.00
chr1	10483	10484	5mC	1000	+	10483	10484	0,0,0	1	100.00
chr1	10488	10489	5mC	1000	+	10488	10489	0,0,0	1	100.00
chr1	10492	10493	5mC	1000	+	10492	10493	0,0,0	1	100.00
chr1	10541	10542	5mC	1000	+	10541	10542	0,0,0	1	100.00
chr1	10562	10563	5mC	1000	+	10562	10563	0,0,0	1	100.00
chr1	10576	10577	5mC	1000	+	10576	10577	0,0,0	1	100.00
chr1	10578	10579	5mC	1000	+	10578	10579	0,0,0	1	0.00
```
#### CNV and PMD

Create folders for output
```
mkdir -p ~/CNV_meth/data/REMORA/PMD_CNV/
```
To merge CNV and PMD methylation results use ~/CNV_meth/Scripts/CNV_meth_REMORA.R and provide as arguments in the following orders:
1) dict_chr.tsv (see above)
2) Remora methylation bed file
3) Segmentation results (hg19 format)
4) output name
5) blacklist bed (hg19 format, not used for this analysis, just leave "NO")
6) whitelist bed (hg38 format, not mandatory): for this analysis we used Partially Methylated Domains from Zhou et al. (2018) 10.1038/s41588-018-0073-4 ~/CNV_meth/Utility/solo_WCGW_inCommonPMDs.hg38.formatted.bed . If you are not using a whitelist leave "NO"
```
/usr/bin/Rscript ~/CNV_meth/Scripts/CNV_meth_REMORA.R ~/CNV_meth/Utility/dict_chr_REMORA.tsv ~/CNV_meth/data/REMORA/sample.remora.hg38.bed  ~/CNV_meth/data/CNV/sample.txt ~/CNV_meth/data/REMORA/PMD_CNV/sample.cnv_meth.R  NO  ~/CNV_meth/Utility/solo_WCGW_inCommonPMDs.hg38.formatted.bed  
```  
#### Average methylation in PMD

Create folders for output
```
mkdir -p ~/CNV_meth/data/REMORA/PMD/
```
To calculate average methylation in PMD use ~/CNV_meth/Scripts/CNV_meth_REMORA.R and provide as arguments in the following orders:
1) dict_chr.tsv (see above)
2) Remora methylation bed file
3) 10MB genome bins (hg19 format), you can find the file used for Katsman et al. in ~/CNV_meth/Utility/hg19_coordinates_total_formatted.txt
4) output name
5) blacklist bed (hg19 format, not used for this analysis, just leave "NO")
6) whitelist bed (hg38 format, not mandatory): for this analysis we used Partially Methylated Domains from Zhou et al. (2018) 10.1038/s41588-018-0073-4 ~/CNV_meth/Utility/solo_WCGW_inCommonPMDs.hg38.formatted.bed . If you are not using a whitelist leave "NO"
```
/usr/bin/Rscript ~/CNV_meth/Scripts/CNV_meth_REMORA.R ~/CNV_meth/Utility/dict_chr_REMORA.tsv ~/CNV_meth/data/REMORA/sample.remora.hg38.bed ~/CNV_meth/Utility/hg19_coordinates_total_formatted.txt ~/CNV_meth/data/REMORA/PMD/sample.cnv_meth.R  NO  ~/CNV_meth/Utility/solo_WCGW_inCommonPMDs.hg38.formatted.bed  
```
#### Genome wide average methylation
Create folders for output
```
mkdir -p ~/CNV_meth/data/REMORA/GENOME/
```

To calculate "genome wide" average methylation use ~/CNV_meth/Scripts/CNV_meth_REMORA.R as before, but use "NO" in the whitelist field
```
/usr/bin/Rscript ~/CNV_meth/Scripts/CNV_meth_REMORA.R ~/CNV_meth/Utility/dict_chr_REMORA.tsv ~/CNV_meth/data/REMORA/sample.remora.hg38.bed  ~/CNV_meth/Utility/hg19_coordinates_total_formatted.txt ~/CNV_meth/data/REMORA/GENOME/sample.cnv_meth.R   NO NO
```
### Deepsignal methylation data

For Deepsignal results perform exactly the same steps but use ~/CNV_meth/Scripts/CNV_meth_DEEPSIGNAL.R and ~/CNV_meth/Utility/dict_chr_DEEPSIGNAL.tsv
put Deepsignal frequency.tsv files in ~/CNV_meth/data/DEEPSIGNAL/ 

```
NC_000005.10	134848203	-	46690055	0.364	0.636	1	0	1	1.0000	CAGCTTCCCGAGTAGCT
NC_000005.10	134848183	-	46690075	0.159	0.841	1	0	1	1.0000	ACTACAGGCGCCAGCCA
NC_000005.10	134848171	-	46690087	0.632	0.368	0	1	1	0.0000	AGCCACCACGCCTGGCT
NC_000001.11	208433305	+	208433305	0.493	0.507	1	0	1	1.0000	AAACAGAGCGTGGAGGA
NC_000018.10	22772251	-	57601033	0.427	0.573	1	0	1	1.0000	TAGTGCTCCGTTGTATA
NC_000010.11	72437128	-	61360293	0.175	0.825	1	0	1	1.0000	GCCCTTAACGCTTTTTC
```
use ~/CNV_meth/data/DEEPSIGNAL/ as output folder


## plots

### CNV

use ~/CNV_meth/scripts/CNV_meth_plot_healthy_and_cancers.R
you can edit the script by specifing the tool used, or directly the folders containing .cnv_meth.R files
```

tool = "DEEPSIGNAL"
tool = "REMORA"

PMD <- makedata(paste("~/CNV_meth/data/", tool ,"/PMD/",sep=""))
GEN <- makedata(paste("~/CNV_meth/data/", tool ,"/GENOME/",sep=""))

```


### 10MB bins

use ~/CNV_meth/Scripts/CNV_meth_plot_cancers.R
you can edit the script by specifing the tool used, or directly the folders containing .cnv_meth.R files
```

tool = "DEEPSIGNAL"
tool = "REMORA"

PMD <- makedata(paste("~/CNV_meth/data/", tool ,"/PMD_CNV/",sep=""))
GEN <- makedata(paste("~/CNV_meth/data/", tool ,"/GENOME/",sep=""))

```


