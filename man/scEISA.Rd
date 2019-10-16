\name{scEISA}
\alias{package scEISA}
\title{Distinguishing the post-transcriptional regulation from transciptional regulation with intron and exon at single cell level
}
\usage{
scEISA(exon_list, intron_list)
}
\description{
The package aims to distinguishing the post-transciptional regulation using a intron matrix and a exon matrix with different treatment.The output will be a csv file in the current directory with gene that are differentially post transcriptionally regulated between 2 treatments.
}
\arguments{
exon_list: A list contains exon reads and treatment.
intron_list: A list contains exon reads and treatment.
}
\examples{
The exon and intron list should include conditions and count matrix:
>exon_list$cont
>['treated', 'treated', 'untreated', 'untreated']

> exon_list$count
> gene_id	bc1	bc2	bc3	bc4
  RPS6	2714	3153	3854394	2513594
  LINC02022	635	7055	2239651	2220470
  RPLP2	1637	2572	1797616	377059
  KCNN1	3981	45172	1688993	1929615
  RPS15	2137	1622	1571169	1002035
  RPL9	5022	4281	1490081	1490487
  ACACA	13022	312404	1301027	2761505
  ST5	14262	215334	1014746	162
  TCF7L1	4431	174011	889380	12340
  AP004609.1	1646	256	888587	0

The output will be:
> scEISA(exon_list, intron_list)
> ENTREZID	logFC	logCPM	LR	PValue	FDR
TSC22D1	TSC22D1	-1.781902171	9.754348291	42.78814926	6.10E-11	7.46E-07
C2orf88	C2orf88	-0.902760176	9.663173855	31.77967988	1.73E-08	9.05E-05
TMPO	TMPO	1.294404205	9.773166734	31.29143891	2.22E-08	9.05E-05
LMNA	LMNA	-0.931711254	9.79633267	16.25656635	5.53E-05	0.135614685
NAGK	NAGK	-0.0003259	9.776056682	16.2524138	5.54E-05	0.135614685
HADHB	HADHB	-0.978333527	9.756083557	14.80829486	0.000119011	0.242603675


}

