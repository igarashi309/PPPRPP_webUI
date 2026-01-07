# Integrated Web Interface for PPR Motif Discovery and Target Prediction in Plants

**[Author Names]**
[Affiliations]

**Corresponding Author:**
Name: [Name]
Address: [Address]
Phone: [Phone]
Email: [Email]

**Running Title:** Unified PPR Motif Analyzer

**Abbreviations:** PPR, pentatricopeptide repeat; GUI, graphical user interface; PLS, P-L-S triplet; ITC, Isothermal Titration Calorimetry.

## Abstract
Pentatricopeptide repeat (PPR) proteins are key regulators of RNA metabolism in plant organelles, involved in essential processes such as RNA editing and splicing. Identifying PPR motifs and predicting their RNA targets remains a complex task requiring multiple computational tools. Here, we present **PPR Motif Analyzer (Unified)**, a standalone software tool that integrates motif discovery and target prediction into a single user-friendly web interface. By combining PPRFinder for accurate motif detection with scoring algorithms based on the PPR code (P-class) and amino acid residues (PLS-class), this tool streamlines the functional analysis of plant PPR proteins. The software provides interactive visualization of motif architecture and RNA binding sites, facilitating rapid hypothesis generation for researchers.

**Key words:** Bioinformatics, Motif discovery, Pentatricopeptide repeat, RNA editing, Web tool.

***

The Pentatricopeptide repeat (PPR) protein family is one of the largest families in land plants, characterized by tandem arrays of a degenerate 35-amino-acid motif. These proteins bind to specific single-stranded RNA sequences in mitochondria and chloroplasts, influencing RNA stability, editing, splicing, and translation (Small and Peeters 2000). Deciphering the specific RNA recognition code of PPR proteins is crucial for understanding organellar gene expression. While previous studies have established the "PPR code" for P-class proteins (Yan et al. 2019) and the specific residue requirements for RNA editing factors in PLS-class proteins (Kobayashi et al. 2019), applying these rules typically requires piecing together disparate command-line tools and custom scripts. This fragmentation poses a technical barrier to many experimental biologists who wish to analyze their proteins of interest.

To address this challenge, we developed **PPR Motif Analyzer (Unified)**, a unified software package that automates the workflow from protein sequence input to target prediction and visualization. The application is built using **Python 3.10+** and utilizes the **Gradio** framework to provide a cross-platform web-based graphical user interface (GUI) that runs on macOS and Linux systems (including WSL for Windows users). The source code is freely available on GitHub.

The software seamlessly integrates motif identification and target scoring. Motif discovery is performed using **PPRFinder** (Gutmann et al. 2020) and **PS_Scan** (ScanProsite) (de Castro et al. 2006). Following motif detection, the tool extracts the specific amino acid residues (typically the 5th and 35th/last positions) that define RNA binding specificity. Two distinct scoring models are implemented to cater to the functional diversity of PPR proteins. For **PLS-class proteins**, which typically function as RNA editing factors, binding scores are calculated using the probability scores for 2-letter amino acid combinations defined by Kobayashi et al. (2019). For **P-class proteins**, binding scores are derived from the dissociation constant ($K_d$) values obtained from Isothermal Titration Calorimetry (ITC) analysis (Yan et al. 2019). In the P-type model, missing $K_d$ values are penalized, and the reciprocal of the $K_d$ ($A = 1/K_d$) is normalized to calculate the base probability ($P_{base}$):
$$P_{base} = \frac{1/K_{d,base}}{\sum_{n \in \{A,C,G,U\}} (1/K_{d,n})}$$

To account for valid structural flexibility in P-type PPR proteins, which may accommodate nucleotide "loop-outs" in the target RNA, the software implements a gap-tolerant scanning algorithm. Users can specify a `max_gap` size (0-5) to investigate alternative alignments. For computational efficiency, the algorithm specifically analyzes site configurations with **no gap** and the **exact requested gap size**, skipping intermediate variations, to identify the optimal binding configuration. Furthermore, to enhance prediction reliability, the final score incorporates a "Motif Confidence Score" ($S_{conf}$) derived from the initial motif search (e.g., HMMER E-values). The final score is a confidence-weighted sum of basis probabilities: $Score_{final} = \sum_{i=1}^{N} (P_{base, i} \times S_{conf, i})$.

The integrated workflow allows users to input a protein FASTA sequence and a target RNA FASTA sequence. The tool automatically executes the pipeline—defining the motif architecture, calculating binding scores across the RNA sequence, and displaying the results as an interactive genome browser track using **Plotly**. Users can manually inspect individual binding sites via a dropdown menu to compare different gap alignment variations (e.g., PLS-type vs P-type alignments), ensuring robustness in site validation. We believe "PPR Motif Analyzer (Unified)" democratizes access to state-of-the-art PPR analysis methods, allowing researchers to generate hypotheses about RNA targets without needing extensive bioinformatics expertise.

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

Gutmann B, Royan S, Schallenberg-Rüdinger M, Lenz H, Castleden IR, McDowell R, Vacher MA, Tonti-Filippini J, Bond CS, Knoop V, Small ID (2020) The Expansion and Diversification of Pentatricopeptide Repeat RNA-Editing Factors in Plants. *Mol Plant* 13: 215–230

Kobayashi T, Yagi Y, Nakamura T (2019) Comprehensive Prediction of Target RNA Editing Sites for PLS-Class PPR Proteins in Arabidopsis thaliana. *Plant Cell Physiol* 60: 862–874

Small ID, Peeters N (2000) The PPR motif−a TPR-related motif prevalent in plant organellar proteins. *Trends Biochem Sci* 25: 46–47

Yan J, Yao Y, Hong S, Yang Y, Shen C, Zhang Q, Zhang D, Zou T, Yin P (2019) Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Res* 47: 3728–3738
