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

  (A) Scripts used for analysing datasets using sci-PDC can be found at <a> <href> https://github.com/reggenlab/Sci-PDC/tree/main/SCI_PDC_CODE </a> </br>


(B) Datasets used in above mentioned scripts can be found at <a> <href> https://github.com/reggenlab/Sci-PDC/tree/main/SCI_PDC_DATA </a> </br>


