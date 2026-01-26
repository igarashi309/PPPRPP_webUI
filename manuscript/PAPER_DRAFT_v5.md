# Integrated Web Interface for PPR Motif Discovery and Target Prediction in Plants

**Keisuke Igarashi**
Graduate School of Agricultural Science, Tohoku University, Sendai, Japan

**Corresponding Author:**
Name: Keisuke Igarashi
Email: [Email]

***

**(Second Page)**

**Running Title:** P-type PPR Prediction Program (PPPRPP)

**Abbreviations:** PPR, pentatricopeptide repeat; GUI, graphical user interface; PLS, P-L-S triplet; ITC, Isothermal Titration Calorimetry.

***

**(Main Text)**

## Abstract
Pentatricopeptide repeat (PPR) proteins are pivotal regulators of organellar RNA metabolism. The PPR family comprises two major subfamilies: the P-class, involved in RNA stabilization and processing, and the PLS-class, primarily functioning as RNA editing factors. While predictive tools for the PLS-class are well-established, software specifically tailored for P-class proteins has been lacking, often forcing researchers to apply suboptimal PLS-type models. I developed **P-type PPR Prediction Program (PPPRPP)**, a cross-platform web application designed for P-type target prediction. PPPRPP integrates a scoring algorithm based on biophysical affinity data ($K_d$) and implements a unique gap-tolerant scanning mechanism to accommodate the flexible binding modes characteristic of P-type proteins. By combining motif discovery (PPRFinder) with these specialized scoring models in a user-friendly GUI, PPPRPP facilitates the functional characterization of the extensive P-class PPR family.

**Key words:** Bioinformatics, Motif discovery, Pentatricopeptide repeat, RNA editing, Web tool.

***

The Pentatricopeptide repeat (PPR) protein family is one of the largest families in land plants, characterized by tandem arrays of a degenerate 35-amino-acid motif. These proteins bind to specific single-stranded RNA sequences in mitochondria and chloroplasts, influencing RNA stability, editing, splicing, and translation (Small and Peeters 2000). Deciphering the specific RNA recognition code of PPR proteins is crucial for understanding organellar gene expression.

Proteins in the PPR family consist of two major subfamilies, defined as P and PLS (Small and Peeters 2000). P-class PPR proteins contain tandem arrays of 35-amino acid PPR motifs. Members of this class have been implicated in RNA stabilization, processing, splicing, and translation. PLS-class proteins contain alternating canonical P-type motifs and variant “long” (L)- and short (S)-type motifs and function mainly in RNA editing. The PLS subfamily is further divided into five subgroups: PLS, E1, E2, E+, and DYW PPR proteins which belong to E, E/E+ and DYW subgroups were almost always involved in RNA editing.

A major challenge in the field has been the discrepancy in tool availability between PPR classes. For PLS-class proteins, which primarily function as RNA editing factors, the "PPR code" and specific residue requirements have been extensively characterized (Kobayashi et al. 2019), and predictive models are widely applied. In contrast, for P-class proteins, which are often involved in RNA stabilization and cleavage, suitable prediction tools have been absent. Although the P-class PPR code has been elucidated through biophysical methods like Isothermal Titration Calorimetry (ITC) (Yan et al. 2019), these experimental affinity values have not been integrated into accessible software. Furthermore, P-type PPR proteins often exhibit a binding mode that tolerates nucleotide "loop-outs" (gaps) in the target RNA, a feature not adequately addressed by existing tools designed for rigid editing factor alignment. This has forced researchers to use PLS-type models for P-type proteins, potentially leading to suboptimal predictions.

The RNA recognition mechanism of PPR proteins has long been proposed as a 1-motif to 1-base correspondence model, but in the P-type PPR subclass, a flexible binding mode where a gap occurs between the number of motifs and the number of recognized bases has been reported (Barkan et al. 2012). This flexibility is particularly evident in studies using maize PPR proteins. Specifically, experimental data from Northern blotting (Pfalz et al. 2009) and gel shift assays (Prikryl et al. 2011) suggest that PPR10 allows a gap of 2-3 bases relative to the total number of PPR motifs in its binding to target RNAs *atpH* and *psaJ*. As an example of a larger gap, CRP1 has been reported to allow a gap of as many as 51 bases, particularly within the 5' untranslated regions of *psaC* and *petA* (Schmitz-Linneweber et al. 2005). The function of CRP1 has been shown to be involved in the protection of transcripts between *petB-petD* and the improvement of translation efficiency of *psaC* (Barkan et al. 1994). These cases suggest that some P-type PPR proteins possess a mechanism to flexibly recognize specific regions of target RNA in addition to the conventional recognition code.

To resolve these issues, I developed **P-type PPR Prediction Program (PPPRPP)**, a comprehensive software package that automates the workflow from protein sequence input to target prediction. The application is built using **Python 3.10+** and utilizes the **Gradio** framework to provide a cross-platform web-based graphical user interface (GUI). The source code is freely available on GitHub (https://github.com/igarashi309/PPPRPP_webUI.git).

The software seamlessly integrates motif identification and target scoring. Motif discovery is performed using **PPRFinder** (Gutmann et al. 2020) and **PS_Scan** (ScanProsite) (de Castro et al. 2006). Motif scanning is a standard method for predicting the positions of known motifs on sequences, and statistical frameworks such as FIMO (Grant et al. 2011) are widely used. These methods estimate the distribution of motif scores at each position based on a background model, calculate corresponding p-values, and evaluate statistical significance by providing q-values through multiple testing correction. However, when performing motif scanning on the entire region of sequences such as plant mitochondrial genomes, which range in length from hundreds of kb to several Mb and contain many repetitive sequences and low-complexity regions, the computational complexity and memory usage increase rapidly, making analysis difficult within the range of realistic computational resources. This problem is caused not only by the increase in score calculation for each position in proportion to the sequence length, but also by the huge number of hits including false positives, making the process of calculating, holding, and correcting p-values and q-values a bottleneck. Therefore, I intentionally did not calculate p-values or q-values, because the statistical assumptions underlying these measures break down when scanning extremely long, repetitive organellar genomes under realistic computational constraints. Instead, I focused on visualizing and comparing the raw motif scores along the genome as an exploratory tool to identify candidate regions for further investigation.

PPPRPP extracts the specific amino acid residues (typically the 5th and 35th/last positions) that define RNA binding specificity and applies one of two user-selectable scoring models. The primary model is the **P-type Model**, based on ITC-derived dissociation constants ($K_d$) for P-class proteins (Yan et al. 2019). The binding scores were derived from the dissociation constant ($K_d$) values obtained from Isothermal Titration Calorimetry (ITC) analysis. In the P-type model, combinations with missing $K_d$ values are penalized by assigning a neutral base probability of **0.25**. The reciprocal of the $K_d$ ($A = 1/K_d$) is normalized to calculate the base probability ($P_{base}$):

$$P_{base} = \frac{1/K_{d,base}}{\sum_{n \in \{A,C,G,U\}} (1/K_{d,n})}$$

The final prediction score can optionally incorporate a "Motif Confidence Score" derived from the motif discovery step (e.g., HMMER bit scores/E-values or PS_Scan scores). If enabled by the user, the base binding probability ($P_{base}$) for each motif is weighted by its confidence score ($S_{conf}$), ensuring that motifs identified with higher certainty contribute more significantly to the target prediction. A secondary **PLS-type Model** is also available, based on probability scores for 2-letter amino acid combinations defined for RNA editing factors (Kobayashi et al. 2019). Binding scores were calculated using the probability scores for 2-letter amino acid combinations (positions 5 and 35/last) defined in the Supplementary Data of Kobayashi et al. (2019). This allows users to compare results or analyze PLS-type proteins within the same interface.

To address the structural flexibility of P-type proteins, the software implements a gap-tolerant scanning algorithm. Users can specify a `max_gap` size (user-definable). The algorithm explicitly evaluates alignments with **no gap** and the **exact requested gap size**, skipping intermediate variations to optimize performance while capturing the characteristic "loop-out" binding mode (e.g., 2-3 nt for PPR10, 51 nt for CRP1).

The tool visualizes results using **Plotly**, providing a **Genome Browser**, an interactive track showing binding scores across the entire RNA molecule; a **Motif Architecture**, a detailed schematic of the protein's motif structure, colored by prediction confidence; and a **Detail Inspector**, a zoomable view of specific binding sites allowing users to inspect the contribution of individual motifs and manually toggle between different gap alignment hypotheses. **Auxiliary Features** include **NCBI Integration** for direct analysis on reference genomes and **IGV Export** for generating standard format files (**BedGraph** for scores, **BED** for binding sites, and **GFF3** for annotations).

To verify the accuracy of the tool, I performed a target prediction for the well-characterized P-type PPR protein **Rf1a** (PPR791) against the *Oryza sativa* mitochondrial genome partial sequence containing the *atp6-orf79* locus. Rf1a is known to promote the cleavage of the *atp6-orf79* cotranscript, silencing the sterility gene *orf79* (Wang et al. 2006). Previous studies have identified the precise binding site of Rf1a within the intergenic region between *atp6* and *orf79* (upstream of *orf79*) (Kazama et al. 2008). Using the default P-type scoring model with `max_gap=0`, PPPRPP successfully accurately reproduced the known binding events (Figure 1). **Figure 1A** shows the Motif Architecture of Rf1a predicted by PPRFinder, displaying the characteristic P-motif repeat structure typical of P-class proteins. **Figure 1B** presents the Genome Browser view, where a distinct, high-scoring peak is clearly visible in the intergenic region between *atp6* and *orf79*, effectively distinguishing the true binding site from background noise. **Figure 1C** provides the Detail Inspector view of this peak, revealing the specific nucleotide-amino acid alignment. Visualizing the amino acids at comparable positions (5 and 35) allows users to verify that the predicted binding site comprises combinations known to have high affinity for the target RNA sequence.

This comprehensive result—from motif detection to specific target site identification—confirms that the software correctly implements the P-class PPR code and accurately predicts experimentally verified targets. The development of PPPRPP addresses a critical need in plant organellar biology: the lack of accessible prediction tools for P-type PPR proteins. By moving beyond the "editing-centric" view and implementing the P-class specific scoring code tailored with $K_d$-based scoring, this tool provides a more physiologically relevant prediction model for proteins involved in RNA stabilization and cleavage. A key innovation of PPPRPP is its explicit handling of "gapped" binding, a phenomenon well-documented in biological studies (e.g., PPR10, CRP1) but largely ignored by rigid alignment tools. While the omission of statistical significance values (p-values) requires users to rely on relative score peaks, this approach is more practical for exploring large, repetitive organellar genomes where traditional motif statistics often break down. I anticipate that PPPRPP will accelerate the functional characterization of the vast number of orphan P-type PPR proteins, shedding light on their diverse roles in post-transcriptional regulation.

***

**Accession numbers**
Sequence data from this article can be found in the GenBank/EMBL data libraries under accession numbers:
*   **AP017386.1**: *Oryza sativa* Indica Group mitochondrial subgenome-2 DNA, complete sequence, strain: BT-CMS.
*   **AB106867.2**: *Oryza sativa* Rf1a mRNA for PPR protein, complete cds.
*   **AAC25599.1**: *Zea mays* CRP1 protein.
*   **X86563**: *Zea mays* chloroplast genome.

**Acknowledgements**
I thank Prof. Kinya Toriyama and laboratory members for their helpful advice and discussions.

**Author Contribution**
K.I. designed the research, performed research, analyzed data, and wrote the paper.

**Funding**
This work was supported by JSPS KAKENHI Grant Number 25K09054.

**Conflict of Interest**
No conflict of interest declared.

**Supplementary Data**
**Supplementary Figure S1.** Gap-tolerant prediction for CRP1.

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

Yan J, Yao Y, Hong S, Yang Y, Shen C, Zhang Q, Zhang D, Zou T, Yin P (2019) Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Res* 47: 3728–3738

Zhelyazkova P, Hammani K, Rojas M, Voelker R, Vargas-Suárez M, Börner T, Barkan A (2012) Protein-mediated protection as the predominant mechanism for defining processed mRNA termini in land plant chloroplasts. *Nucleic Acids Res* 40: 3092-3105

**Figures**

**Figure 1. Application of PPPRPP to the P-type protein Rf1a (PPR791).**
(A) Motif architecture of Rf1a predicted by PPRFinder. The protein consists of 18 P-class PPR motifs (colored blocks).
(B) Genome Browser view of the *Oryza sativa* mitochondrial genome region (133kb–134kb) containing the *atp6* and *orf79* genes. The P-type scoring model (max_gap=0) identifies a distinct high-scoring peak (Score > 180) in the intergenic region, corresponding to the known Rf1a binding site.
(C) Detail Inspector view of the peak binding site. The alignment shows the 1-to-1 correspondence between the 5th/Last amino acid residues of the PPR motifs (top rows) and the target RNA sequence (bottom row), consistent with the standard P-class recognition code.

**Supplementary Figure S1. Gap-tolerant prediction for CRP1.**
(A) Motif architecture of the P-type protein CRP1.
(B) Genome Browser view of the *Zea mays* chloroplast genome near the *psaC* locus. Analysis using the P-type model with `max_gap=51` reveals a significant binding peak.
(C) Detail Inspector view of the predicted CRP1 binding site. The alignment visualizes the "loop-out" mode, where a 51-nucleotide gap (highlighted) is accommodated between specific motifs to achieve high-affinity binding, consistent with the experimentally reported stabilization mechanism of the *psaC* 5' UTR.
