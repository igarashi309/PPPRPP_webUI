# Integrated Web Interface for PPR Motif Discovery and Target Prediction in Plants

## Abstract
Pentatricopeptide repeat (PPR) proteins are key regulators of RNA metadata in plant organelles, involved in essential processes such as RNA editing and splicing. Identifying PPR motifs and predicting their RNA targets remains a complex task requiring multiple computational tools. Here, we present **PPR Motif Analyzer (Unified)**, a standalone software tool that integrates motif discovery and target prediction into a single user-friendly web interface. By combining `PPRFinder` for accurate motif detection with scoring algorithms based on the PPR code (P-class) and amino acid residues (PLS-class), this tool streamlines the functional analysis of plant PPR proteins. The software provides interactive visualization of motif architecture and RNA binding sites, facilitating rapid hypothesis generation for researchers.

## Introduction
The Pentatricopeptide repeat (PPR) protein family is one of the largest families in land plants, characterized by tandem arrays of a degenerate 35-amino-acid motif. These proteins bind to specific single-stranded RNA sequences in mitochondria and chloroplasts, influencing RNA stability, editing, splicing, and translation. Deciphering the binding code of PPR proteins is crucial for understanding organellar gene expression.

While previous studies have established the "PPR code" for P-class proteins (Yan et al., 2019) and the specific residue requirements for RNA editing factors in PLS-class proteins (Kobayashi et al., 2019), applying these rules typically requires piecing together disparate command-line tools and scripts. This fragmentation poses a barrier to experimental biologists.

To address this, we developed a unified software package that automates the workflow from protein sequence input to target prediction visualization. 

## Materials and Methods

### Software Architecture
The application is built using **Python 3.10+** and utilizes the **Gradio** framework for its web-based graphical user interface (GUI). It runs on macOS and Linux systems (including WSL for Windows users).

### Motif Discovery
Motif identification was performed using **PPRFinder** (Gutmann et al., 2020) and **PS_Scan** (ScanProsite) (de Castro et al., 2006).

### Target Prediction Algorithms
The software implements two distinct scoring models tailored to different functional classes of PPR proteins:
1.  **PLS-type (Kobayashi et al., 2019)**: Optimized for PLS-class proteins involved in RNA editing. Binding scores were calculated using the probability scores for 2-letter amino acid combinations (positions 5 and 35/last) defined in the Supplementary Data of Kobayashi et al. (2019).
2.  **P-type (Yan et al., 2019)**: Designed for P-class proteins. Binding scores were derived from the dissociation constant ($K_d$) values obtained from Isothermal Titration Calorimetry (ITC) analysis by Yan et al. (2019). The binding probability ($P_{base}$) for each base was calculated as follows:
    1.  Missing $K_d$ values were assigned a penalty value of 99,999 $\mu M$.
    2.  The reciprocal of the $K_d$ value was taken to represent affinity ($A = 1/K_d$).
    3.  Values were normalized to sum to 1 across the four bases:
        $$P_{base} = \frac{1/K_{d,base}}{\sum_{n \in \{A,C,G,U\}} (1/K_{d,n})}$$

### Visualization
Binding results are visualized using **Plotly**, enabling interactive exploration of the genome. Users can view binding scores across the entire RNA sequence or zoom into specific genomic ranges. The tool also generates BedGraph and BED files for compatibility with the **Integrative Genomics Viewer (IGV)**.

## Results and Features

### Unified Workflow
Users provide a protein FASTA sequence and a target RNA FASTA sequence. The tool performs the following steps automatically:
1.  executes `hmmsearch` and `PPRFinder` to define the motif architecture.
2.  Extracts the "PPR code" (5th and last amino acids) from each motif.
3.  Scans the RNA sequence to calculate binding scores for every position using the selected scoring model (PLS-type or P-type).
4.  Displays the results as an interactive genome browser track.

### Interactive Interface
The GUI allows users to:
-   Toggle between "PPRFinder" and "PS_Scan" for motif discovery.
-   Select the appropriate scoring method ("PLS-type" or "P-type").
-   Inspect individual binding sites to see the contribution of each motif to the total score.
-   Download results in CSV format for further analysis.
-   Directly fetch reference genomes and annotations from NCBI using Accession IDs.

## Conclusion
PPR Motif Analyzer (Unified) provides a streamlined, accessible solution for researchers studying plant RNA metabolism. By encapsulating complex computational steps into a single interface, it democratizes access to state-of-the-art PPR analysis methods.

## Availability
The source code is available on GitHub: https://github.com/igarashi309/PPPRPP_webUI

## References

1.  **Kobayashi, T., Yagi, Y., & Nakamura, T. (2019).** Comprehensive Prediction of Target RNA Editing Sites for PLS-Class PPR Proteins in Arabidopsis thaliana. *Plant and Cell Physiology*, 60(4), 862–874. [https://doi.org/10.1093/pcp/pcy251](https://doi.org/10.1093/pcp/pcy251)
2.  **Yan, J., Yao, Y., Hong, S., Yang, Y., Shen, C., Zhang, Q., Zhang, D., Zou, T., & Yin, P. (2019).** Delineation of pentatricopeptide repeat codes for target RNA prediction. *Nucleic Acids Research*, 47(7), 3728–3738. [https://doi.org/10.1093/nar/gkz075](https://doi.org/10.1093/nar/gkz075)
3.  **Gutmann, B., Royan, S., Schallenberg-Rüdinger, M., Lenz, H., Castleden, I. R., McDowell, R., Vacher, M. A., Tonti-Filippini, J., Bond, C. S., Knoop, V., & Small, I. D. (2020).** The Expansion and Diversification of Pentatricopeptide Repeat RNA-Editing Factors in Plants. *Molecular Plant*, 13(2), 215–230. [https://doi.org/10.1016/j.molp.2019.11.002](https://doi.org/10.1016/j.molp.2019.11.002)
4.  **de Castro, E., Sigrist, C. J. A., Gattiker, A., Bulliard, V., Langendijk-Genevaux, P. S., Gasteiger, E., Bairoch, A., & Hulo, N. (2006).** ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. *Nucleic Acids Research*, 34(suppl_2), W362–W365. [https://doi.org/10.1093/nar/gkl124](https://doi.org/10.1093/nar/gkl124)
