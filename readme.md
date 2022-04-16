<h1> Single-Cell based Inference of association between Pathway, Disease, and Cell-type (sci-PDC) </h1>

sci-PDC is a method for performing single-cell expression-based inference of relationships between pathway, disease, and cell type. This is the main repository for sci-PDC.

<h2>Basic workflow </h2>
<ul>
(I) Gene-set enrichment analysis of single-cell expression profiles. </br>
(II) Co-occurrence analysis for estimating:
  <li>Correlation (Spearman) </li>
 <li>Specificity score (based on pre-defined null model)</li> 
 <li>Rank specificity (based on rank normalisation by cell atlases)</li> 
(III) Data filteration </br>
(IV) Learning network using Prababilistic Graph Model (PGM):
  <li>Bayesian Network (Directed Graph)</li>
  <li>Markov Network (Undirected Graph)</li> 
(V)PubMed based validation
  </ul>
  
  
<h4> The repository contains :-- </br> </h4>
(A) Scripts in  <a> <href> https://github.com/reggenlab/Sci-PDC/tree/main/SCI_PDC_CODE </a>  which corresponds to: </br>

Code for single-cell expression profile based enrichment analysis <a> <href> https://github.com/reggenlab/Sci-PDC/blob/main/SCI_PDC_CODE/enrichment_score.R </a> </br> 
</br>
Three main case studies (using diabetes mellitus as disease example): </br>
(I) Analysisng human pancreatic beta cells for disease-pathway associations <a> <href>https://github.com/reggenlab/Sci-PDC/blob/main/SCI_PDC_CODE/case_study1.R </a> </br>
(II) Exploring variability across species in terms of disease-pathway associations in human pancreatic beta cells <a> <href>https://github.com/reggenlab/Sci-PDC/blob/main/SCI_PDC_CODE/case_study2.R </a> </br>
(III) Investigating effect of age in terms of disease-pathway associations for human pancreatic beta cells <a> <href> https://github.com/reggenlab/Sci-PDC/blob/main/SCI_PDC_CODE/case_study3.R </a> </br>
</br>
Comparison of different approaches for disease-pathway association analysis <a> <href> https://github.com/reggenlab/Sci-PDC/blob/main/SCI_PDC_CODE/comparasion_plot.R </a> </br>

(B) Datasets used in above mentioned scripts can be found at <a> <href> https://github.com/reggenlab/Sci-PDC/tree/main/SCI_PDC_DATA </a> </br>


