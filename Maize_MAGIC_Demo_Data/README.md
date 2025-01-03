Test commonly used functions for **Maize** MAGIC data from [MaizeMAGIC](https://github.com/rqtl/qtl2data/tree/main/MaizeMAGIC) (for ***R/qtl2***).  An associated publication is [Dell'Acqua et al. 2015](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0716-z).  Relative to [ArabMAGIC](https://github.com/rqtl/qtl2data/tree/main/ArabMAGIC), there are additional variables for *maize_magic_crossinfo.csv* and *maize_magic_phenocovar.csv*.

`run_R_qtl2_demo-Maize-HeightOnly.R` - similar to [run_R_qtl2_demo-Height.R](https://github.com/cwarden45/BPSC234_Term_Paper-QTL/blob/main/Arabidopsis_MAGIC_Demo_Data/run_R_qtl2_demo-Height.R) for *Arabidopsis*, but use Maize MAGIC data (and "PH" for height).

`run_R_qtl2_demo-Maize-HeightPhenocovar.R` - similar to `run_R_qtl2_demo-Maize-HeightOnly.R`, but provide full traits and only return "PH" results.

## *R/qtl2* results (for *PH* or "*P*lant-*H*eight")

![](MaizeMAGIC_Demo-HeightOnly-Combined_LOD.png)

For the Haley-Knott Regression, the **95% permutation threshold** for the LOD Score was **7.4** (after setting the random seed to 0).

That can be seen in the separate plot below (where **7.4** is noticeably higher than a LOD of 3 or 4):

![](MaizeMAGIC_Demo-HeightOnly-HK_LOD.png)

Similar results can be shown for height (*PH*) when the **other traits are considered as covariates**:

<table>
  <tbody>
    <tr>
      <th align="center">Phenotype</th>
      <th align="center">Permutation Threshold</th>
    </tr>
    <tr>
	<td align="left"><b>PS</b></br>(<i>"pollen shed"</i>)</td>
	<td align="center">7.44</td>
    </tr>
    <tr>
	<td align="left"><b>PH</b></br>(<i>"plant height"</i>)</td>
	<td align="center">7.4</td>
    </tr>
    <tr>
	<td align="left"><b>EH</b></br>(<i>"ear height"</i>)</td>
	<td align="center">7.72</td>
    </tr>
    <tr>
	<td align="left"><b>GYrad</b></br>(<i>"transformed grain yield"</i>)</td>
	<td align="center">7.49</td>
    </tr>
</tbody>
</table>

![](MaizeMAGIC_Demo-HeightPhenocovar-Combined_LOD.png)

![](MaizeMAGIC_Demo-HeightPhenocovar-HK_LOD.png)


## Testing *R/qtl2* `est_map()` Function

If I try to run `test_R_qtl2_demo-Maize-est_map.R`, then I encounter the following error message in the result object:

***est_map not yet implemented for general RIL.***
