\name{edgeR_RQLF_withCont}
\alias{edgeR_RQLF_withCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential expression analysis using edgeR_RQLF with a contrast matrix
}
\description{
The function aims to explore the differential expression analysis based on edgeR_RQLF along with a contrast matrix creating the interaction term of study of interest.To identify the post transcriptional regulation from transcriptional changes, the interaction term is delta exon - delta intron.
}
\usage{
edgeR_RQLF_withCont(exon-and-intron-List)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exon-and-intron-List}{A list contains both intron and exon reads and treatments and a contrast matrix which creates interaction term for study of interest.
}}

{
\examples{
The exon and intron list should include conditions:
>exon-and-intron-List$cont
>['ex_cont1', 'ex_cont2', 'in_cont1', 'in_cont2']

> exon-and-intron-List$count
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

The contrast matrix should be like:
> L$cont.matrix
          Contrasts
Levels     in_cont1_cont2 ex_cont1_cont2 exonMin_intron
  ex_cont1              0              1              1
  ex_cont2              0             -1             -1
  in_cont1              1              0             -1
  in_cont2             -1              0              1

The output will be:
> edgeR_RQLF_withCont(exon-and-intron-List, contrast)
> ENTREZID	logFC	logCPM	LR	PValue	FDR
TSC22D1	TSC22D1	-1.781902171	9.754348291	42.78814926	6.10E-11	7.46E-07
C2orf88	C2orf88	-0.902760176	9.663173855	31.77967988	1.73E-08	9.05E-05
TMPO	TMPO	1.294404205	9.773166734	31.29143891	2.22E-08	9.05E-05
LMNA	LMNA	-0.931711254	9.79633267	16.25656635	5.53E-05	0.135614685
NAGK	NAGK	-0.0003259	9.776056682	16.2524138	5.54E-05	0.135614685
HADHB	HADHB	-0.978333527	9.756083557	14.80829486	0.000119011	0.242603675


}



