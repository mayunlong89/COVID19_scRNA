# COVID19_scRNA
Here we submitted codes and methods on the manuscript entitled "Single cell sequencing and GWAS analyses reveal genetics-influenced CD16+monocytes and memory CD8+T cells involved in severe COVID-19".

# Abstract
Background: Understanding the host genetic architecture and viral immunity contribute to develop effective vaccines and therapeutics for controlling the COVID-19 pandemic. Alterations of immune responses in peripheral blood mononuclear cells play a crucial role in the detrimental progression of COVID-19. To date, the effects of host genetic factors on peripheral immune response for severe COVID-19 remain largely unknown. 
Methods: We collected the GWAS summary statistics on hospitalized COVID-19 with large-scale samples (N = 969,689) from the COVID-19 Host Genetics Consortium, and four scRNA-seq datasets with large-scale cells (N = 606,534) from the GEO database. We constructed a powerful computational framework to characterize the host genetics-influenced immune cell subpopulations by integrating GWAS summary statistics with scRNA-seq data.
Results: Total 34 risk genes and 10 biological pathways were significantly associated with severe COVID-19 based on meta-GWAS summary statistics. The number of highly-expressed genetics-risk genes and pathways were remarkably elevated with increased severity of COVID-19. Integrative analysis of GWAS and single cell sequencing identified the three cell subpopulations of CD16+monocytes, megakaryocytes, and memory CD8+T cells were significantly enriched by severe COVID-19-related genetic components. Notably, three causal risk genes of CCR1, CXCR6, and ABO were specifically expressed in these three cell types, respectively. CCR1+CD16+monocytes and ABO+ megakaryocytes with significant up-regulated genes of S100A12, S100A8, S100A9, and IFITM1 confer risk to the cytokine storms among severe patients. CXCR6+ memory CD8+ T cells exhibit a notable polyfunctionality of multiple immunologic features, including elevation of proliferation, migration, and chemotaxis. Additionally, we observed a prominently increased number of cell-cell interactions of CCR1+ CD16+monocytes and CXCR6+ memory CD8+T cells in severe patients compared to normal controls among both PBMCs and lung tissues, and elevated interactions of epithelial cells contribute to enhance the resident to lung airway for against COVID-19 infection.
Conclusions: We uncover a major genetics-modulated immunological shift between mild and severe infection, including an increase in up-regulated genetic-risk genes, excessive secreted inflammatory cytokines, and functional immune cell subsets contributing high risk to severity. Our study provides novel insights in parsing the host genetics-influenced immune cells for severe COVID-19.

# Introduction
Accumulating evidence have suggested alterations of immune responses in peripheral blood mononuclear cells (PBMCs) play a crucial role in the detrimental progression of COVID-19. A growing number of GWASs have identified numerous significant genetic variants associated with COVID-19 susceptibility and severity. Many earlier GWASs have shown that complex genetic dysregulations of peripheral immune cells with highly selective effects on the risk of immune-related diseases at the subcellular level. However, the effect of these genetic determinants on the peripheral immune cells for severe COVID-19 remains largely unknown. In light of there is no comprehensive study for revealing the genetically regulatory effects of peripheral immune cells on severe COVID-19, the present study is the first integrative genomic analysis by combining genetic information from GWAS with scRNA-seq data to genetically pinpoint immune cell types implicated in the etiology of severe COVID-19. 

# Methods
## 1. Single cell RNA-seq data on severe COVID-19  
In the current study, we downloaded four independent scRNA-seq datasets on COVID-19 and its severity in PBMC and BALF from the ArrayExpress database (Dataset #1: the accession number is E-MTAB-9357) from Su et al. study [11], and the Gene Expression Omnibus (GEO) database (Dataset #2: the accession number is GSE149689 from Lee et al. study [20], Dataset #3: the accession number is GSE150861 from Guo et al. study [12], and Dataset #4: the accession number is GSE158055 [9]). For dataset #1, this dataset contained 270 peripheral blood samples including 254 samples with different COVID-19 severity (i.e., mild N = 109, moderate N = 102, and severe N = 50) and 16 healthy controls for scRNA-seq analysis. For the dataset #2, there were eight patients with COVID-19 of varying clinical severity, including asymptomatic, mild, and severe, and four healthy controls with PBMCs. As for the dataset #3, there were five peripheral blood samples from two severe COVID-19 patients at three different time points during tocilizumab treatment, containing two different stages: severe stage and remission stage. With regard to the dataset #4, there were 12 BALF samples including three moderate and nine severe patients collected from lung tissues. For all datasets, the sample collection process underwent Institutional Review Board review and approval at the institutions where samples were originally collected. The COVID-19 severity was qualified by using the World Health Organization (WHO) ordinal scale (WOS), the National Early Warning Score (NEWS), or the Diagnosis and Treatment of COVID-19 (Trail Version 6). Single-cell transcriptomes for these four datasets were gathered by using the 10× Genomics scRNA-seq platform. 

## 2. GWAS summary statistics from the COVID-19 Host Genetic Consortium
  The meta-GWAS summary data on severe COVID-19 round 4 (B2_ALL, Susceptibility [Hospitalized COVID-19 vs. Population]) were downloaded from the official website of the COVID-19 Host Genetic Consortium [23] (https://www.covid19hg.org/; analyzed file named: “COVID19_HGI_B2_ALL_leave_23andme_20201020.txt.gz”; released date of October 4 2020). There were 7,885 hospitalized COVID-19 patients and 961,804 control participants from 21 independent contributing studies. The vast majority of participants in these contributing studies were of European ancestry (93%). The meta-GWAS summary statistics contained P values, Wald statistic, inverse-variance meta-analyzed log Odds Ratio (OR) and related standard errors. The 1,000 Genomes Project European Phase 3 [37] were used as a panel for pruning. Results from 23&Me cohort GWAS summary statistics were excluded from our current analysis. By filtering genetic variants without RefSNP number in the Human Genome reference builds 37, there were 9,368,170 genetic variants included with a major allele frequency (MAF) threshold of 0.0001 and the imputation score filter of 0.6. We used the qqman R package for figuring the Manhattan plot to visualize the meta-GWAS analysis results. The web-based software of LocusZoom [38] was utilized to visualize the regional association plots for identified risk loci (http://locuszoom.sph.umich.edu/).
  
 ## 3. Scripts:
  In the present sutyd, we leveraged numerous bioinformatics tools: linux-based tools incluidng MAGMA, S-MultiXcan, R-based tools including Rolypoly and permutation, and web-access tools inclduing the WEB-based Gene SeT AnaLysis Toolkit (WebGestalt; http://www.webgestalt.org) [42], the PhenoScanner V2 (http://www.phenoscanner.medschl.cam.ac.uk/) [45],  the Open Target Genetics (OTG, https://genetics.opentargets.org/) [46], STRING(v11.0, https://string-db.org/)[51], STITCH (v5.0, http://stitch.embl.de/)[53],ChEMBL (v2.6, https://www.ebi.ac.uk/chembl/) [54], and DGIdb database (https://www.dgidb.org/druggable_gene_categories). 
  In order to ensure our peers could follow our analyses, we have deposited the codes and methods in the current github, as the following example:
```
#compute rolypoly
######################
library("rolypoly")
library("data.table")
#
index<-c("normal","mild","moderate","severe")
lapply(index,function(x){
  file_n<-paste0("/share/pub/dengcy/Singlecell/COVID19/data/Rploy_",x,"_cell.txt")
  merge_scexpr<-read.delim(file_n,sep = " ")
  colnames(merge_scexpr)<-annotation$V2
  merge_scexpr<-merge_scexpr[apply(merge_scexpr,1,sum)!=0,]
  #create the annotation files
  gene_name<-intersect(rownames(merge_scexpr),geneid_df1$label)
  geneid_df1<-geneid_df1[geneid_df1$label %in% gene_name,]
  merge_scexpr<-merge_scexpr[gene_name,]
  geneid_df1<-geneid_df1[!duplicated(geneid_df1$label),]
#############################################
  file_na<-paste0("roly_",x,"_pre.RData")
  save(geneid_df1,merge_scexpr,file=file_na)
  })

 ld_path <- "/share/pub/dengcy/Singlecell/COVID19/data/LD"

#sim_block_annotation$label<-rownames(merge_scexprc2)[1:1000]
#Rploy_remission_GSE.txt
 rolypoly_result <- rolypoly_roll(
   gwas_data = COVID19_GWAS_autosomes_maf2,
   block_annotation = geneid_df1,
   block_data = merge_scexpr,
   ld_folder =ld_path,
   bootstrap_iters = 100
  )
  save(rolypoly_result,file = "/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/rolypoly_mild_cell.RData")
```
# Reference
1.	Dong E, Du H, Gardner L: An interactive web-based dashboard to track COVID-19 in real time. Lancet Infect Dis 2020, 20:533-534.
2.	Wu Z, McGoogan JM: Characteristics of and Important Lessons From the Coronavirus Disease 2019 (COVID-19) Outbreak in China: Summary of a Report of 72 314 Cases From the Chinese Center for Disease Control and Prevention. JAMA 2020.
3.	Berlin DA, Gulick RM, Martinez FJ: Severe Covid-19. N Engl J Med 2020.
4.	Richardson S, Hirsch JS, Narasimhan M, Crawford JM, McGinn T, Davidson KW, Barnaby DP, Becker LB, Chelico JD, Cohen SL, et al: Presenting Characteristics, Comorbidities, and Outcomes Among 5700 Patients Hospitalized With COVID-19 in the New York City Area. JAMA 2020, 323:2052-2059.
5.	Guan WJ, Ni ZY, Hu Y, Liang WH, Ou CQ, He JX, Liu L, Shan H, Lei CL, Hui DSC, et al: Clinical Characteristics of Coronavirus Disease 2019 in China. N Engl J Med 2020, 382:1708-1720.
6.	Xu L, Ma Y, Yuan J, Zhang Y, Wang H, Zhang G, Tu C, Lu X, Li J, Xiong Y, et al: COVID-19 Quarantine Reveals Behavioral Changes Effect on Myopia Progression. Ophthalmology 2021.
7.	Pedersen SF, Ho YC: SARS-CoV-2: a storm is raging. J Clin Invest 2020, 130:2202-2205.
8.	Takahashi T, Ellingson MK, Wong P, Israelow B, Lucas C, Klein J, Silva J, Mao T, Oh JE, Tokuyama M, et al: Sex differences in immune responses that underlie COVID-19 disease outcomes. Nature 2020, 588:315-320.
9.	Chen G, Wu D, Guo W, Cao Y, Huang D, Wang H, Wang T, Zhang X, Chen H, Yu H, et al: Clinical and immunological features of severe and moderate coronavirus disease 2019. J Clin Invest 2020, 130:2620-2629.
10.	Su Y, Chen D, Yuan D, Lausted C, Choi J, Dai CL, Voillet V, Duvvuri VR, Scherler K, Troisch P, et al: Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell 2020, 183:1479-1495.e1420.
11.	Guo C, Li B, Ma H, Wang X, Cai P, Yu Q, Zhu L, Jin L, Jiang C, Fang J, et al: Single-cell analysis of two severe COVID-19 patients reveals a monocyte-associated and tocilizumab-responding cytokine storm. Nat Commun 2020, 11:3924.
12.	Ren X, Wen W, Fan X, Hou W, Su B, Cai P, Li J, Liu Y, Tang F, Zhang F, et al: COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas. Cell 2021.
13.	Wen W, Su W, Tang H, Le W, Zhang X, Zheng Y, Liu X, Xie L, Li J, Ye J, et al: Immune cell profiling of COVID-19 patients in the recovery stage by single-cell sequencing. Cell Discov 2020, 6:31.
14.	Zhang JY, Wang XM, Xing X, Xu Z, Zhang C, Song JW, Fan X, Xia P, Fu JL, Wang SY, et al: Single-cell landscape of immunological responses in patients with COVID-19. Nat Immunol 2020, 21:1107-1118.
15.	Chua RL, Lukassen S, Trump S, Hennig BP, Wendisch D, Pott F, Debnath O, Thürmann L, Kurth F, Völker MT, et al: COVID-19 severity correlates with airway epithelium-immune cell interactions identified by single-cell analysis. Nat Biotechnol 2020, 38:970-979.
16.	Delorey TM, Ziegler CGK, Heimberg G, Normand R, Yang Y, Segerstolpe A, Abbondanza D, Fleming SJ, Subramanian A, Montoro DT, et al: A single-cell and spatial atlas of autopsy tissues reveals pathology and cellular targets of SARS-CoV-2. bioRxiv 2021.
17.	Silvin A, Chapuis N, Dunsmore G, Goubet AG, Dubuisson A, Derosa L, Almire C, Hénon C, Kosmider O, Droin N, et al: Elevated Calprotectin and Abnormal Myeloid Cell Subsets Discriminate Severe from Mild COVID-19. Cell 2020, 182:1401-1418 e1418.
18.	Schulte-Schrepping J, Reusch N, Paclik D, Baßler K, Schlickeiser S, Zhang B, Krämer B, Krammer T, Brumhard S, Bonaguro L, et al: Severe COVID-19 Is Marked by a Dysregulated Myeloid Cell Compartment. Cell 2020, 182:1419-1440 e1423.
19.	Lee JS, Park S, Jeong HW, Ahn JY, Choi SJ, Lee H, Choi B, Nam SK, Sa M, Kwon JS, et al: Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19. Sci Immunol 2020, 5.
20.	Cao X: COVID-19: immunopathology and its implications for therapy. Nat Rev Immunol 2020, 20:269-270.
21.	Del Valle DM, Kim-Schulze S, Huang HH, Beckmann ND, Nirenberg S, Wang B, Lavin Y, Swartz TH, Madduri D, Stock A, et al: An inflammatory cytokine signature predicts COVID-19 severity and survival. Nat Med 2020, 26:1636-1643.
22.	Arunachalam PS, Wimmers F, Mok CKP, Perera R, Scott M, Hagan T, Sigal N, Feng Y, Bristow L, Tak-Yin Tsang O, et al: Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans. Science 2020, 369:1210-1220.
23.	The COVID-19 Host Genetics Initiative, a global initiative to elucidate the role of host genetic factors in susceptibility and severity of the SARS-CoV-2 virus pandemic. Eur J Hum Genet 2020, 28:715-718.
24.	Ellinghaus D, Degenhardt F, Bujanda L, Buti M, Albillos A, Invernizzi P, Fernández J, Prati D, Baselli G, Asselta R, et al: Genomewide Association Study of Severe Covid-19 with Respiratory Failure. N Engl J Med 2020.
25.	Pathak GA, Singh K, Miller-Fleming TW, Wendt FR, Ehsan N, Hou K, Johnson R, Lu Z, Gopalan S, Yengo L, et al: Integrative analyses identify susceptibility genes underlying COVID-19 hospitalization. medRxiv 2020:2020.2012.2007.20245308.
26.	Shelton JF, Shastri AJ, Ye C, Weldon CH, Filshtein-Somnez T, Coker D, Symons A, Esparza-Gordillo J, Aslibekyan S, Auton A: Trans-ethnic analysis reveals genetic and non-genetic associations with COVID-19 susceptibility and severity. medRxiv 2020:2020.2009.2004.20188318.
27.	Roberts GHL, Park DS, Coignet MV, McCurdy SR, Knight SC, Partha R, Rhead B, Zhang M, Berkowitz N, Haug Baltzell AK, et al: AncestryDNA COVID-19 Host Genetic Study Identifies Three Novel Loci. medRxiv 2020:2020.2010.2006.20205864.
28.	Zhou S, Butler-Laporte G, Nakanishi T, Morrison DR, Afilalo J, Afilalo M, Laurent L, Pietzner M, Kerrison N, Zhao K, et al: A Neanderthal OAS1 isoform protects individuals of European ancestry against COVID-19 susceptibility and severity. Nat Med 2021, 27:659-667.
29.	Pairo-Castineira E, Clohisey S, Klaric L, Bretherick AD, Rawlik K, Pasko D, Walker S, Parkinson N, Fourman MH, Russell CD, et al: Genetic mechanisms of critical illness in COVID-19. Nature 2021, 591:92-98.
30.	Ma Y, Huang Y, Zhao S, Yao Y, Zhang Y, Qu J, Wu N, Su J: Integrative Genomics Analysis Reveals a 21q22.11 Locus Contributing Risk to COVID-19. Hum Mol Genet 2021.
31.	Gaziano L, Giambartolomei C, Pereira AC, Gaulton A, Posner DC, Swanson SA, Ho YL, Iyengar SK, Kosik NM, Vujkovic M, et al: Actionable druggable genome-wide Mendelian randomization identifies repurposing opportunities for COVID-19. Nat Med 2021, 27:668-676.
32.	Orrù V, Steri M, Sidore C, Marongiu M, Serra V, Olla S, Sole G, Lai S, Dei M, Mulas A, et al: Complex genetic signatures in immune cells underlie autoimmunity and inform therapy. Nat Genet 2020, 52:1036-1045.
33.	Roederer M, Quaye L, Mangino M, Beddall MH, Mahnke Y, Chattopadhyay P, Tosi I, Napolitano L, Terranova Barberio M, Menni C, et al: The genetic architecture of the human immune system: a bioresource for autoimmunity and disease pathogenesis. Cell 2015, 161:387-403.
34.	Patin E, Hasan M, Bergstedt J, Rouilly V, Libri V, Urrutia A, Alanio C, Scepanovic P, Hammer C, Jönsson F, et al: Natural variation in the parameters of innate immune cells is preferentially driven by genetic factors. Nat Immunol 2018, 19:302-314.
35.	Aguirre-Gamboa R, Joosten I, Urbano PCM, van der Molen RG, van Rijssen E, van Cranenbroek B, Oosting M, Smeekens S, Jaeger M, Zorro M, et al: Differential Effects of Environmental and Genetic Factors on T and B Cell Immune Traits. Cell Rep 2016, 17:2474-2487.
36.	Organization WH: COVID-19 Theapeutic Trial Synopsis. (World Health Organization) 2020.
37.	Auton A, Brooks LD, Durbin RM, Garrison EP, Kang HM, Korbel JO, Marchini JL, McCarthy S, McVean GA, Abecasis GR: A global reference for human genetic variation. Nature 2015, 526:68-74.
38.	Pruim RJ, Welch RP, Sanna S, Teslovich TM, Chines PS, Gliedt TP, Boehnke M, Abecasis GR, Willer CJ: LocusZoom: regional visualization of genome-wide association scan results. Bioinformatics 2010, 26:2336-2337.
39.	Calderon D, Bhaskar A, Knowles DA, Golan D, Raj T, Fu AQ, Pritchard JK: Inferring Relevant Cell Types for Complex Traits by Using Single-Cell Gene Expression. Am J Hum Genet 2017, 101:686-699.
40.	Chang CC, Chow CC, Tellier LC, Vattikuti S, Purcell SM, Lee JJ: Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience 2015, 4:7.
41.	de Leeuw CA, Mooij JM, Heskes T, Posthuma D: MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput Biol 2015, 11:e1004219.
42.	Wang J, Duncan D, Shi Z, Zhang B: WEB-based GEne SeT AnaLysis Toolkit (WebGestalt): update 2013. Nucleic Acids Res 2013, 41:W77-83.
43.	Kanehisa M, Goto S: KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res 2000, 28:27-30.
44.	Ma Y, Li J, Xu Y, Wang Y, Yao Y, Liu Q, Wang M, Zhao X, Fan R, Chen J, et al: Identification of 34 genes conferring genetic and pharmacological risk for the comorbidity of schizophrenia and smoking behaviors. Aging (Albany NY) 2020, 12:2169-2225.
45.	Staley JR, Blackshaw J, Kamat MA, Ellis S, Surendran P, Sun BB, Paul DS, Freitag D, Burgess S, Danesh J, et al: PhenoScanner: a database of human genotype-phenotype associations. Bioinformatics 2016, 32:3207-3209.
46.	Ghoussaini M, Mountjoy E, Carmona M, Peat G, Schmidt EM, Hercules A, Fumis L, Miranda A, Carvalho-Silva D, Buniello A, et al: Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics. Nucleic Acids Res 2021, 49:D1311-d1320.
47.	Barbeira AN, Dickinson SP, Bonazzola R, Zheng J, Wheeler HE, Torres JM, Torstenson ES, Shah KP, Garcia T, Edwards TL, et al: Exploring the phenotypic consequences of tissue specific gene expression variation inferred from GWAS summary statistics. Nat Commun 2018, 9:1825.
48.	Barbeira AN, Pividori M, Zheng J, Wheeler HE, Nicolae DL, Im HK: Integrating predicted transcriptome from multiple tissues improves association detection. PLoS Genet 2019, 15:e1007889.
49.	Ma X, Wang P, Xu G, Yu F, Ma Y: Integrative genomics analysis of various omics data and networks identify risk genes and variants vulnerable to childhood-onset asthma. BMC Med Genomics 2020, 13:123.
50.	Xu M, Li J, Xiao Z, Lou J, Pan X, Ma Y: Integrative genomics analysis identifies promising SNPs and genes implicated in tuberculosis risk based on multiple omics datasets. Aging (Albany NY) 2020, 12:19173-19220.
51.	Ma Y, Huang Y, Zhao S, Yao Y, Zhang Y, Qu J, Wu N, Su J: Integrative Genomics Analysis Reveals a Novel 21q22.11 Locus Contributing to Susceptibility of COVID-19. medRxiv 2020:2020.2009.2016.20195685.
52.	von Mering C, Huynen M, Jaeggi D, Schmidt S, Bork P, Snel B: STRING: a database of predicted functional associations between proteins. Nucleic Acids Res 2003, 31:258-261.
53.	Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T: Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res 2003, 13:2498-2504.
54.	Szklarczyk D, Santos A, von Mering C, Jensen LJ, Bork P, Kuhn M: STITCH 5: augmenting protein-chemical interaction networks with tissue and affinity data. Nucleic Acids Res 2016, 44:D380-384.
55.	Gaulton A, Hersey A, Nowotka M, Bento AP, Chambers J, Mendez D, Mutowo P, Atkinson F, Bellis LJ, Cibrián-Uhalte E, et al: The ChEMBL database in 2017. Nucleic Acids Res 2017, 45:D945-d954.
56.	Mathew D, Giles JR, Baxter AE, Oldridge DA, Greenplate AR, Wu JE, Alanio C, Kuri-Cervantes L, Pampena MB, D'Andrea K, et al: Deep immune profiling of COVID-19 patients reveals distinct immunotypes with therapeutic implications. Science 2020, 369.
57.	Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, Mauck WM, 3rd, Hao Y, Stoeckius M, Smibert P, Satija R: Comprehensive Integration of Single-Cell Data. Cell 2019, 177:1888-1902.e1821.
58.	Pairo-Castineira E, Clohisey S, Klaric L, Bretherick AD, Rawlik K, Pasko D, Walker S, Parkinson N, Fourman MH, Russell CD, et al: Genetic mechanisms of critical illness in Covid-19. Nature 2020.
59.	Battle A, Brown CD, Engelhardt BE, Montgomery SB: Genetic effects on gene expression across human tissues. Nature 2017, 550:204-213.
60.	Wang Q, Chen R, Cheng F, Wei Q, Ji Y, Yang H, Zhong X, Tao R, Wen Z, Sutcliffe JS, et al: A Bayesian framework that integrates multi-omics data and gene networks predicts risk genes from schizophrenia GWAS data. Nat Neurosci 2019, 22:691-699.
61.	Ma Y, Li MD: Establishment of a Strong Link Between Smoking and Cancer Pathogenesis through DNA Methylation Analysis. Sci Rep 2017, 7:1811.
62.	Zhu Z, Zhang F, Hu H, Bakshi A, Robinson MR, Powell JE, Montgomery GW, Goddard ME, Wray NR, Visscher PM, Yang J: Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets. Nat Genet 2016, 48:481-487.
63.	He X, Fuller CK, Song Y, Meng Q, Zhang B, Yang X, Li H: Sherlock: detecting gene-disease associations by matching patterns of expression QTL and GWAS. Am J Hum Genet 2013, 92:667-680.
64.	Zhou P, Yang XL, Wang XG, Hu B, Zhang L, Zhang W, Si HR, Zhu Y, Li B, Huang CL, et al: A pneumonia outbreak associated with a new coronavirus of probable bat origin. Nature 2020, 579:270-273.
65.	Lu R, Zhao X, Li J, Niu P, Yang B, Wu H, Wang W, Song H, Huang B, Zhu N, et al: Genomic characterisation and epidemiology of 2019 novel coronavirus: implications for virus origins and receptor binding. Lancet 2020, 395:565-574.
66.	Wilk AJ, Rustagi A, Zhao NQ, Roque J, Martínez-Colón GJ, McKechnie JL, Ivison GT, Ranganath T, Vergara R, Hollis T, et al: A single-cell atlas of the peripheral immune response in patients with severe COVID-19. Nat Med 2020, 26:1070-1076.
67.	Hadjadj J, Yatim N, Barnabei L, Corneau A, Boussier J, Smith N, Péré H, Charbit B, Bondet V, Chenevier-Gobeaux C, et al: Impaired type I interferon activity and inflammatory responses in severe COVID-19 patients. Science 2020, 369:718-724.
68.	Merad M, Martin JC: Pathological inflammation in patients with COVID-19: a key role for monocytes and macrophages. Nat Rev Immunol 2020, 20:355-362.
69.	Dimopoulos G, de Mast Q, Markou N, Theodorakopoulou M, Komnos A, Mouktaroudi M, Netea MG, Spyridopoulos T, Verheggen RJ, Hoogerwerf J, et al: Favorable Anakinra Responses in Severe Covid-19 Patients with Secondary Hemophagocytic Lymphohistiocytosis. Cell Host Microbe 2020, 28:117-123.e111.
70.	Jamilloux Y, Henry T, Belot A, Viel S, Fauter M, El Jammal T, Walzer T, François B, Sève P: Should we stimulate or suppress immune responses in COVID-19? Cytokine and anti-cytokine interventions. Autoimmun Rev 2020, 19:102567.

