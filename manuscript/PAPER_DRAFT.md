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
Pentatricopeptide repeat (PPR) proteins are key regulators of RNA metabolism in plant organelles, involved in essential processes such as RNA editing and splicing. While predictive models for PLS-class PPR proteins (RNA editing factors) are well-established, tools specifically designed for P-class PPR proteins (involved in RNA stabilization and processing) have been lacking. Consequently, researchers have often been forced to apply PLS-class scoring models to P-class proteins, despite functional and structural differences. To address this gap, we developed **P-type PPR Prediction Program (PPPRPP)**, a unified software package tailored for P-type PPR target prediction. The program implements a scoring algorithm based on binding affinity ($K_d$) data from Isothermal Titration Calorimetry (ITC) of P-class proteins and incorporates a gap-tolerant scanning mechanism to account for the structural flexibility observed in P-type RNA recognition. The tool also supports PLS-class scoring, enabling comparative analysis. By integrating motif discovery (PPRFinder) and these specialized scoring models into a user-friendly web interface, PPPRPP facilitates the specific and effective identification of P-type PPR targets.

**Key words:** Bioinformatics, Motif discovery, Pentatricopeptide repeat, RNA editing, Web tool.

***

The Pentatricopeptide repeat (PPR) protein family is one of the largest families in land plants, characterized by tandem arrays of a degenerate 35-amino-acid motif. These proteins bind to specific single-stranded RNA sequences in mitochondria and chloroplasts, influencing RNA stability, editing, splicing, and translation (Small and Peeters 2000). Deciphering the specific RNA recognition code of PPR proteins is crucial for understanding organellar gene expression.

A major challenge in the field has been the discrepancy in tool availability between PPR classes. For PLS-class proteins, which primarily function as RNA editing factors, the "PPR code" and specific residue requirements have been extensively characterized (Kobayashi et al. 2019), and predictive models are widely applied. In contrast, for P-class proteins, which are often involved in RNA stabilization and cleavage, suitable prediction tools have been absent. Although the P-class PPR code has been elucidated through biophysical methods like Isothermal Titration Calorimetry (ITC) (Yan et al. 2019), these experimental affinity values have not been integrated into accessible software. Furthermore, P-type PPR proteins often exhibit a binding mode that tolerates nucleotide "loop-outs" (gaps) in the target RNA, a feature not adequately addressed by existing tools designed for the rigid alignment of editing factors. This has forced researchers to use PLS-type models for P-type proteins, potentially leading to suboptimal predictions.

To resolve these issues, we developed **P-type PPR Prediction Program (PPPRPP)**, a comprehensive software package that automates the workflow from protein sequence input to target prediction. The application is built using **Python 3.10+** and utilizes the **Gradio** framework to provide a cross-platform web-based graphical user interface (GUI). The source code is freely available on GitHub.

The software seamlessly integrates motif identification and target scoring. Motif discovery is performed using **PPRFinder** (Gutmann et al. 2020) and **PS_Scan** (ScanProsite) (de Castro et al. 2006). Statistical frameworks such as FIMO (Grant et al. 2011) are widely used for motif scanning but face significant hurdles when applied to plant organellar genomes. These genomes are large (hundreds of kb to Mb) and rich in repetitive sequences, causing computational complexity and false discovery rates to explode. Calculating and correcting p-values/q-values under these conditions creates a bottleneck and raises issues with statistical consistency across fragmented analyses. Therefore, we intentionally omitted p-value/q-value calculations. Instead, PPPRPP focuses on visualizing raw motif scores across the genome as an exploratory tool to identify candidate regions based on relative affinity peaks.

**Scoring Algorithms**
PPPRPP extracts the specific amino acid residues (typically the 5th and 35th/last positions) that define RNA binding specificity and applies one of two user-selectable scoring models:
1.  **P-type Model (Primary):** Based on ITC-derived dissociation constants ($K_d$) for P-class proteins (Yan et al. 2019). The score represents the normalized binding probability ($P_{base} \propto 1/K_d$). Missing $K_d$ values are penalized.
2.  **PLS-type Model (Secondary):** Based on probability scores for 2-letter amino acid combinations defined for RNA editing factors (Kobayashi et al. 2019). This allows users to compare results or analyze PLS-type proteins within the same interface.

**Gap-Tolerant Scanning**
To address the structural flexibility of P-type proteins, the software implements a gap-tolerant scanning algorithm. Users can specify a `max_gap` size (0-5). The algorithm explicitly evaluates alignments with **no gap** and the **exact requested gap size**, skipping intermediate variations to optimize performance while capturing the characteristic "loop-out" binding mode.

**Visualization and Interface**
The tool visualizes results using **Plotly**, providing:
-   **Genome Browser:** An interactive track showing binding scores across the entire RNA molecule.
-   **Motif Architecture:** A detailed schematic of the protein's motif structure, colored by prediction confidence, with labels for motif types (P, L, S) and their specific indices (e.g., #1, #2).
-   **Detail Inspector:** A zoomable view of specific binding sites, allowing users to inspect the contribution of individual motifs and manually toggle between different gap alignment hypotheses.

**Validation**
To verify the accuracy of the tool, we performed a target prediction for the well-characterized P-type PPR protein **Rf1a** (PPR791) against the *Oryza sativa* mitochondrial genome partial sequence containing the *atp6-orf79* locus. Rf1a is known to bind to the *atp6-orf79* cotranscript to promote cleavage (Wang et al. 2006). Using the default P-type scoring model with `max_gap=0`, PPPRPP successfully identified the known binding site upstream of *atp6* as a prominent peak (Figure 1). This result confirms that the software correctly implements the P-class PPR code and scoring logic.

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
de Castro E, Sigrist CJA, Gattiker A, Bulliard V, Langendijk-Genevaux PS, Gasteiger E, Bairoch A, Hulo N (2006) ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. *Nucleic Acids Res* 34: W362–W365

Grant CE, Bailey TL, Noble WS (2011) FIMO: scanning for occurrences of a given motif. *Bioinformatics* 27: 1017–1018

Gutmann B, Royan S, Schallenberg-Rüdinger M, Lenz H, Castleden IR, McDowell R, Vacher MA, Tonti-Filippini J, Bond CS, Knoop V, Small ID (2020) The Expansion and Diversification of Pentatricopeptide Repeat RNA-Editing Factors in Plants. *Mol Plant* 13: 215–230

Kobayashi T, Yagi Y, Nakamura T (2019) Comprehensive Prediction of Target RNA Editing Sites for PLS-Class PPR Proteins in Arabidopsis thaliana. *Plant Cell Physiol* 60: 862–874

Small ID, Peeters N (2000) The PPR motif−a TPR-related motif prevalent in plant organellar proteins. *Trends Biochem Sci* 25: 46–47

Wang Z, Zou Y, Li X, Zhang Q, Chen L, Wu H, Su D, Guo Y, Guo D, Long Y, Ge Y, Liu R (2006) Cytoplasmic male sterility of rice with boro II cytoplasm is associated with a cytotoxic peptide and is restored by two related PPR motif genes via distinct modes of mRNA silencing. *Plant Cell* 18: 676-687

Yan J, Yao Y, Hong S, Yang Y, Shen C, Zhang Q, Zhang D, Zou T, Yin P (2019) Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Res* 47: 3728–3738
