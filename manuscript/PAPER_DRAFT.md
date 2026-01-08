# Integrated Web Interface for PPR Motif Discovery and Target Prediction in Plants

**[Author Names]**
[Affiliations]

**Corresponding Author:**
Name: [Name]
Address: [Address]
Phone: [Phone]
Email: [Email]

**Running Title:** P-type PPR Prediction Program (PPPRPP)

**Abbreviations:** PPR, pentatricopeptide repeat; GUI, graphical user interface; PLS, P-L-S triplet; ITC, Isothermal Titration Calorimetry.

## Abstract
Pentatricopeptide repeat (PPR) proteins are key regulators of RNA metabolism in plant organelles, involved in essential processes such as RNA editing and splicing. While predictive models for PLS-class PPR proteins (RNA editing factors) are well-established, tools specifically designed for P-class PPR proteins (involved in RNA stabilization and processing) have been lacking. Consequently, researchers have often been forced to apply PLS-class scoring models to P-class proteins, despite functional and structural differences. To address this gap, we developed **P-type PPR Prediction Program (PPPRPP)**, a unified software package tailored for P-type PPR target prediction. The program implements a scoring algorithm based on binding affinity ($K_d$) data from Isothermal Titration Calorimetry (ITC) of P-class proteins and incorporates a gap-tolerant scanning mechanism to account for the structural flexibility observed in P-type RNA recognition. The tool also supports PLS-class scoring, enabling comparative analysis. By integrating motif discovery (PPRFinder) and these specialized scoring models into a user-friendly web interface, PPPRPP supports the exploration of potential P-type PPR targets.

**Key words:** Bioinformatics, Motif discovery, Pentatricopeptide repeat, RNA editing, Web tool.

***

The Pentatricopeptide repeat (PPR) protein family is one of the largest families in land plants, characterized by tandem arrays of a degenerate 35-amino-acid motif. These proteins bind to specific single-stranded RNA sequences in mitochondria and chloroplasts, influencing RNA stability, editing, splicing, and translation (Small and Peeters 2000). Deciphering the specific RNA recognition code of PPR proteins is crucial for understanding organellar gene expression.

A major challenge in the field has been the discrepancy in tool availability between PPR classes. For PLS-class proteins, which primarily function as RNA editing factors, the "PPR code" and specific residue requirements have been extensively characterized (Kobayashi et al. 2019), and predictive models are widely applied. In contrast, for P-class proteins, which are often involved in RNA stabilization and cleavage, suitable prediction tools have been absent. Although the P-class PPR code has been elucidated through biophysical methods like Isothermal Titration Calorimetry (ITC) (Yan et al. 2019), these experimental affinity values have not been integrated into accessible software. Furthermore, P-type PPR proteins often exhibit a binding mode that tolerates nucleotide "loop-outs" (gaps) in the target RNA, a feature not adequately addressed by existing tools designed for the rigid alignment of editing factors. This has forced researchers to use PLS-type models for P-type proteins, potentially leading to suboptimal predictions.

The RNA recognition mechanism of PPR proteins has long been proposed as a 1-motif to 1-base correspondence model, but in the P-type PPR subclass, a flexible binding mode where a gap occurs between the number of motifs and the number of recognized bases has been reported (Barkan et al. 2012). This flexibility is particularly evident in studies using maize PPR proteins. Specifically, experimental data from Northern blotting (Pfalz et al. 2009) and gel shift assays (Prikryl et al. 2011) suggest that PPR10 allows a gap of 2-3 bases relative to the total number of PPR motifs in its binding to target RNAs *atpH* and *psaJ*. As an example of a larger gap, CRP1 has been reported to allow a gap of as many as 9 bases upon binding to target RNA. The function of CRP1 has been shown to be involved in the protection of transcripts between *petB-petD* (Barkan et al. 1994; Zhelyazkova et al. 2012) and the improvement of translation efficiency of *psaC* (Barkan et al. 1994; Schmitz-Linneweber et al. 2005), and functional analysis of these supports the existence of binding that does not strictly correspond to the number of motifs. These cases suggest that some P-type PPR proteins possess a mechanism to flexibly recognize specific regions of target RNA in addition to the conventional recognition code.

To resolve these issues, we developed **P-type PPR Prediction Program (PPPRPP)**, a comprehensive software package that automates the workflow from protein sequence input to target prediction. The application is built using **Python 3.10+** and utilizes the **Gradio** framework to provide a cross-platform web-based graphical user interface (GUI). The source code is freely available on GitHub.

The software seamlessly integrates motif identification and target scoring. Motif discovery is performed using **PPRFinder** (Gutmann et al. 2020) and **PS_Scan** (ScanProsite) (de Castro et al. 2006). Motif scanning is a standard method for predicting the positions of known motifs on sequences, and statistical frameworks such as FIMO (Grant et al. 2011) are widely used. These methods estimate the distribution of motif scores at each position based on a background model, calculate corresponding p-values, and evaluate statistical significance by providing q-values through multiple testing correction. However, when performing motif scanning on the entire region of sequences such as plant mitochondrial genomes, which range in length from hundreds of kb to several Mb and contain many repetitive sequences and low-complexity regions, the computational complexity and memory usage increase rapidly, making analysis difficult within the range of realistic computational resources. This problem is caused not only by the increase in score calculation for each position in proportion to the sequence length, but also by the huge number of hits including false positives, making the process of calculating, holding, and correcting p-values and q-values a bottleneck. Furthermore, if the analysis target is divided to reduce the computational load, p-values and q-values are calculated based on different test populations for each analysis unit, which raises the problem that direct comparison between runs is not statistically valid. Therefore, there is an inherent trade-off between "statistically consistent motif prediction for the entire huge sequence" and "realistically feasible computational resources".

For these reasons, in the "search for motif-like sequences targeting the entire plant mitochondrial genome" which is the purpose of this study, it is difficult to apply the existing statistical framework as it is, and it is necessary to explicitly deal with the gap between computational resource constraints and statistical consistency. Consequently, we intentionally did not calculate p-values or q-values, because the statistical assumptions underlying these measures break down when scanning extremely long, repetitive organellar genomes under realistic computational constraints. Instead, we focused on visualizing and comparing the raw motif scores along the genome as an exploratory tool to identify candidate regions for further investigation.

**Scoring Algorithms**
PPPRPP extracts the specific amino acid residues (typically the 5th and 35th/last positions) that define RNA binding specificity and applies one of two user-selectable scoring models35: 1.  **P-type Model (Primary):** Based on ITC-derived dissociation constants ($K_d$) for P-class proteins (Yan et al. 2019). The binding scores were derived from the dissociation constant ($K_d$) values obtained from Isothermal Titration Calorimetry (ITC) analysis. In the P-type model, combinations with missing $K_d$ values are penalized by assigning a neutral base probability of **0.25**. The reciprocal of the $K_d$ ($A = 1/K_d$) is normalized to calculate the base probability ($P_{base}$):
36: $$P_{base} = \frac{1/K_{d,base}}{\sum_{n \in \{A,C,G,U\}} (1/K_{d,n})}$$
37: The final prediction score can optionally incorporate a "Motif Confidence Score" derived from the motif discovery step (e.g., HMMER bit scores/E-values or PS_Scan scores). If enabled by the user, the base binding probability ($P_{base}$) for each motif is weighted by its confidence score ($S_{conf}$), ensuring that motifs identified with higher certainty contribute more significantly to the target prediction:
$$Score_{final} = \sum_{i=1}^{N} (P_{base, i} \times S_{conf, i})$$

2.  **PLS-type Model (Secondary):** Based on probability scores for 2-letter amino acid combinations defined for RNA editing factors (Kobayashi et al. 2019). Binding scores were calculated using the probability scores for 2-letter amino acid combinations (positions 5 and 35/last) defined in the Supplementary Data of Kobayashi et al. (2019). This allows users to compare results or analyze PLS-type proteins within the same interface.

42: **Gap-Tolerant Scanning**
43: To address the structural flexibility of P-type proteins, the software implements a gap-tolerant scanning algorithm. Users can specify a `max_gap` size (user-definable, up to 20 nucleotides). The algorithm explicitly evaluates alignments with **no gap** and the **exact requested gap size**, skipping intermediate variations to optimize performance while capturing the characteristic "loop-out" binding mode (e.g., 2-3 nt for PPR10, 9 nt for CRP1).

**Visualization and Interface**
The tool visualizes results using **Plotly**, providing:
-   **Genome Browser:** An interactive track showing binding scores across the entire RNA molecule.
-   **Motif Architecture:** A detailed schematic of the protein's motif structure, colored by prediction confidence, with labels for motif types (P, L, S) and their specific indices (e.g., #1, #2).
49: -   **Detail Inspector:** A zoomable view of specific binding sites, allowing users to inspect the contribution of individual motifs and manually toggle between different gap alignment hypotheses.
50: 
51: **Auxiliary Features**
52: -   **NCBI Integration:** Users can perform analysis directly on reference genomes by retrieving sequences via NCBI Accession IDs.
53: -   **IGV Export:** The tool generates standard format files (**BedGraph** for scores, **BED** for binding sites, and **GFF3** for annotations) that can be exported for visualization in external genome browsers like the Integrative Genomics Viewer (IGV).

**Validation**
To verify the accuracy of the tool, we performed a target prediction for the well-characterized P-type PPR protein **Rf1a** (PPR791) against the *Oryza sativa* mitochondrial genome partial sequence containing the *atp6-orf79* locus. Rf1a is known to promote the cleavage of the *atp6-orf79* cotranscript, silencing the sterility gene *orf79*. Previous studies have identified the precise binding site of Rf1a within the intergenic region between *atp6* and *orf79* (upstream of *orf79*) (Kazama et al. 2008). Using the default P-type scoring model with `max_gap=0`, PPPRPP successfully identified a prominent peak that exactly matches the reported Rf1a binding site (Figure 1). This result confirms that the software correctly implements the P-class PPR code and accurately predicts experimentally verified targets.

**Accession numbers**
Sequence data from this article can be found in the GenBank/EMBL data libraries under accession numbers [Accession Numbers].

**Acknowledgements**
We thank [Name] for [Contribution].

**Author Contribution**
[Initials] designed the research; [Initials] performed research; [Initials] analyzed data; and [Initials] wrote the paper.

**Funding**
This work was supported by [Funding Source] (Grant No. [Number] to [Initials]).

**Conflict of Interest**
No conflict of interest declared.

**Supplementary Data**
**Supplementary Figure S1.** [Title of Figure S1]
**Supplementary Table S1.** [Title of Table S1]

**References**
Barkan A, Rojas M, Fujii S, Yap A, Chong YS, Bond CS, Small I (2012) A combinatorial amino acid code for RNA recognition by pentatricopeptide repeat proteins. *Proc Natl Acad Sci USA* 109: 15535-15540

Barkan A, Walker M, Nolasco M, Johnson D (1994) A nuclear mutation in maize blocks the processing and translation of several chloroplast mRNAs and provides evidence for the differential translation of alternative mRNA forms. *EMBO J* 13: 3170-3181

de Castro E, Sigrist CJA, Gattiker A, Bulliard V, Langendijk-Genevaux PS, Gasteiger E, Bairoch A, Hulo N (2006) ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. *Nucleic Acids Res* 34: W362–W365

Grant CE, Bailey TL, Noble WS (2011) FIMO: scanning for occurrences of a given motif. *Bioinformatics* 27: 1017–1018

Gutmann B, Royan S, Schallenberg-Rüdinger M, Lenz H, Castleden IR, McDowell R, Vacher MA, Tonti-Filippini J, Bond CS, Knoop V, Small ID (2020) The Expansion and Diversification of Pentatricopeptide Repeat RNA-Editing Factors in Plants. *Mol Plant* 13: 215–230

Kobayashi T, Yagi Y, Nakamura T (2019) Comprehensive Prediction of Target RNA Editing Sites for PLS-Class PPR Proteins in Arabidopsis thaliana. *Plant Cell Physiol* 60: 862–874

Kazama T, Nakamura T, Watanabe M, Sugita M, Toriyama K (2008) Suppression mechanism of mitochondrial ORF79 accumulation by Rf1 protein in BT-type cytoplasmic male sterile rice. *Plant J* 55: 619-628

Pfalz J, Bayraktar OA, Prikryl J, Barkan A (2009) Site-specific binding of a PPR protein defines and stabilizes 5' and 3' mRNA termini in chloroplasts. *EMBO J* 28: 2042-2052

Prikryl J, Rojas M, Schuster G, Barkan A (2011) Mechanism of RNA stabilization and translational activation by a pentatricopeptide repeat protein. *Proc Natl Acad Sci USA* 108: 415-420

Schmitz-Linneweber C, Williams-Carrier R, Barkan A (2005) RNA immunoprecipitation and microarray analysis show a chloroplast pentatricopeptide repeat protein to be associated with the 5′ region of mRNAs whose translation it activates. *Plant Cell* 17: 2791-2804

Small ID, Peeters N (2000) The PPR motif−a TPR-related motif prevalent in plant organellar proteins. *Trends Biochem Sci* 25: 46–47

Wang Z, Zou Y, Li X, Zhang Q, Chen L, Wu H, Su D, Guo Y, Guo D, Long Y, Ge Y, Liu R (2006) Cytoplasmic male sterility of rice with boro II cytoplasm is associated with a cytotoxic peptide and is restored by two related PPR motif genes via distinct modes of mRNA silencing. *Plant Cell* 18: 676-687

Zhelyazkova P, Hammani K, Rojas M, Voelker R, Vargas-Suárez M, Börner T, Barkan A (2012) Protein-mediated protection as the predominant mechanism for defining processed mRNA termini in land plant chloroplasts. *Nucleic Acids Res* 40: 3092-3105

Yan J, Yao Y, Hong S, Yang Y, Shen C, Zhang Q, Zhang D, Zou T, Yin P (2019) Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Res* 47: 3728–3738
