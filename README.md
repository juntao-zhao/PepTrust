# PepTrust: A Probabilistic Framework for High-Confidence Peptide Identification

Liquid chromatography-tandem mass spectrometry (LC-MS/MS)-based peptide identification is fundamental to proteomics. When evaluating peptide identification results, the lack of ground truth in real biological samples is a significant challenge. This hampers downstream analyses, such as protein identification and quantification, as well as follow-up experiment design, including mechanistic and drug target studies in biomedical research. To address this challenge, we propose a method named PepTrust to estimate confidence scores of peptide identification results without requiring ground truth, enhancing reliability in complex proteomics datasets.

PepTrust is based on a truth discovery framework [1], modeling peptide identification methods (e.g., Comet, X! Tandem, MS-GF+) as sources that provide peptide identity claims for spectra. It integrates these claims through a probabilistic approach, estimating the trustworthiness of each source based on consistency with respect to outputs of multiple combination methods. By evaluating the likelihood of observed spectra given the peptide identification results, PepTrust assigns confidence scores to peptide identifications, ranking them to identify those with the highest reliability. This approach capitalizes on consensus across different combination methods while adjusting for their varying reliability, offering a robust alternative to traditional ground-truth-dependent validation approaches.

Preliminary evaluations on public proteomics datasets demonstrate that PepTrust effectively isolates high-confidence peptides. On a sub-pool of the PRIDE dataset PDX004732 with known ground truth [2], PepTrust achieved 70.29% precision, surpassing the 64.91% average across individual search methods. This framework strengthens analyses in ground-truth-absent scenarios, facilitating more reliable downstream applications.

[1] Fanget al. ACMTIST 11.6 (2020): 1-24.

[2] Zolg, et al. Nature Methods 14.3 (2017): 259-262.

## Key Words
Mass Spectrometry, Peptide identification, ground-truth-absent evaluation, probabilistic framework, truth discovery
# PepTrust
