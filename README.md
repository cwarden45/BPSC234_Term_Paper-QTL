The *[Bloom et al. 2019](https://elifesciences.org/articles/49212)* *eLife* publication / dataset was selected due the combination of [processed genotype data](https://www.dropbox.com/scl/fo/o764opdwcly2cpo061sqc/ACQBKlJPPzo2VS4Marwcr2o?rlkey=y01ub37eysvj0kb6df1dj2swh&e=1&dl=0), [public code](https://github.com/joshsbloom/yeast-16-parents), [raw data deposit](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA549760/), and an [earlier preprint discussion](http://disq.us/p/216tpbo).

`stats_to_select_trait.R` - script used to guide selection of **trait** to select for re-analysis for project (based upoon [elife-49212-fig4-data1-v2.xlsx](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDkyMTIvZWxpZmUtNDkyMTItZmlnNC1kYXRhMS12Mi54bHM-/elife-49212-fig4-data1-v2.xls?_hash=XMfZAGDNF6CLi7La1FnFaGthQcgTKD52MgylCT9sPjc%3D) for ***[Figure 4b](https://iiif.elifesciences.org/lax/49212%2Felife-49212-fig4-v2.tif/full/1500,/0/default.jpg)***) as well as [NIHMS544073-supplement-01.xls](https://docs.google.com/spreadsheets/d/1LLmAJKEacVwooooJxWB-rliZSMS1ONeh/edit?usp=drive_link&ouid=113994953300724355569&rtpof=true&sd=true) (which is within *RData-20241012T200813Z-001.zip*, if the full set of data is downloaded for [RData](https://drive.google.com/drive/folders/1nIV-DgJQ0qPnyqj1UWbN1FkpfBlnJH1k)).
 - Decide to focus on **lactate** trait.


`stats_to_select_cross.R` - script used to guide selection of **cross** (or **crosses**) to select for re-analysis for project (based upoon [elife-49212-fig2-data1-v2.xlsx](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDkyMTIvZWxpZmUtNDkyMTItZmlnMy1kYXRhMi12Mi54bHM-/elife-49212-fig3-data2-v2.xls?_hash=YDQCI094FNfj1wfMNr1qda86EUlLSrrS4cPrhGyIBLk%3D) for ***[Figure 3a](https://iiif.elifesciences.org/lax/49212%2Felife-49212-fig3-v2.tif/full/1500,/0/default.jpg)***) as well as [elife-49212-fig4-data1-v2.xlsx](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDkyMTIvZWxpZmUtNDkyMTItZmlnNC1kYXRhMS12Mi54bHM-/elife-49212-fig4-data1-v2.xls?_hash=XMfZAGDNF6CLi7La1FnFaGthQcgTKD52MgylCT9sPjc%3D) (with was also used for `stats_to_select_trait.R`).
 - **MKT1** is the gene that is associate with the *lactate* trait in this study as well as *"High temperature growth, Sporulation efficiency, Gene expression, DNA damage, petite formation, ethanol tolerance, dipropyldopamine resistance"* traits in muliple earlier studies (*"Steinmetz et al. 2002, Deutschbauer and Davis 2005, Sinha et al. 2006, Smith and Kruglyak 2006, Demogines et al. 2008, Dimitrov et al. 2009, Swinnen et al. 2011, Kim and Fay 2009"*).
 - **MKT1** is identified to contain a causal variant related to the *lactate* trait with a variant inherited from **BYa**.
 - Based considering variability affecting the trait within a cross and crosses that involve *BYa*, I thought that I should consider **Cross "A"** (***BYxRM***)
 - In particular, the *G30D* (**chrXIV_467219_A_G**) variant in *MKT1* is defined as a *"Rare"* variant in this study.
 - In this study, the variant is associated with **many other traits** when inherited from *BYa* (*Caffeine (15mM)*, *EtOH*, *Glycerol*, *Trehalose*, *Magnesium_Chloride*, *Raffinose*, *YPD (37)*, *Sucrose*, *Xylose*, *Copper_Sulfate*, *Galactose*, *Maltose*, and *Congo_Red (75ug/mL)*): the trend is usually a decrease, but it is an increase for *Caffeine*.  The variant is **also** associated with *Mannose* when the shared parent is *RMx* (for the decreasing trend).  This may be a reason to also consider the **Cross "375"** (***M22xBY***), especially if there is some focus on this variant.

The yeast segregants are **haploid**, so that may have some influence on the methods that can be tested (in addition to the experimental design).

`run_R_qtl2-Lactate-Model1.R` - script to reformat input files and run ***R/qtl2*** analysis based upon [user guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html) example (but with a `haploid` cross design).

***Counts for chrXIV_467219_A_G Genotypes***

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center">Geno1</th>
      <th align="center">Geno2</th>
    </tr>
    <tr>
	<td align="left">Cross "A" (<i>BYxRM</i>)</td>
	<td align="center">468</td>
	<td align="center">483</td>
    </tr>
    <tr>
	<td align="left">Cross "375" (<i>M22xBY</i>)</td>
	<td align="center"></td>
	<td align="center"></td>
    </tr>
</tbody>
</table>

**Preiliminary *R/qtl2* Lactate Results**:

||Cross "A" (*BYxRM*)|Cross "375" (*M22xBY*)|
| --- | --- |--- |
|Model 1 |![](A_BYxRM_Lactate.Model1_LODall.png)|![](A_M22xBY_Lactate.Model1_LODall.png)|

<table>
  <tbody>
    <tr>
      <th align="center"></th>
      <th align="center"><b><i>R/qtl2</i></b><br>(haploid,<br>no covariates)</th>
    </tr>
    <tr>
	<td align="left">Cross "A" (<i>BYxRM</i>)</td>
	<td align="center"></td>
    </tr>
    <tr>
	<td colspan=2"></td>
    </tr>
    <tr>
	<td align="left">Cross "375" (<i>M22xBY</i>)</td>
	<td align="center"></td>
    </tr>
</tbody>
</table>

Even though there is a trial version, I did **not** look further into the applicable experimental designs because there is a license for long-term use of *[MapQTL](https://www.kyazma.nl/index.php/MapQTL/)*.

Other publications / datasets used for general preparatation (for the term project and course content) include *[Gandolfi et al. 2018](https://www.nature.com/articles/s41598-018-25438-0)* (which has been partially used for earlier analysis of domestic cat Whole Genome Sequencing data [here](https://github.com/cwarden45/Bastu_Cat_Genome) and [here](https://github.com/cwarden45/MissYvonne_Cat_Genome), and there is more specific discussion of the genetic map creation/calculation in *[Li et al. 2016](https://academic.oup.com/g3journal/article/6/6/1607/6029961)*), *[Yengo et al. 2022](https://www.nature.com/articles/s41586-022-05275-y)* (saturation of common variants for human height), and *[Fei_et_al_2022](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03470-7)* (with maize SLAF-Seq [raw data deposit](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA778629)).  For example, even if not used for re-analysis of this specfiic dataset, *[Kontou and Bagos 2024](https://biodatamining.biomedcentral.com/articles/10.1186/s13040-024-00385-x)* provides a review (with some emphasis on analysis that considered summary statistics), and *[Hofmeister et al. 2023](https://www.nature.com/articles/s41588-023-01415-w)* is an example of a relatively recent update / development for a genetic analysis method (*SHAPEIT5*).  In terms of my individiual general research interests, I also thought *[Herniter et al. 2019](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.01346/full)* was a good fit.  Based upon this early preparation, I thought content for the ***R/qtl2*** (*[Broman et al. 2019](https://academic.oup.com/genetics/article/211/2/495/5931111)*) could be relevant for this paper/project, and that is why I started with testing that package for analysis.  *Slides with "Notes" created during the preparation process have been uploaded as PDF files.*
