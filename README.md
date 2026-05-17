# A reproduction of the chronology building procedure in Buntgen et al. (2020)

In Buntgen et al. (2020), the authors use a sparse set of ensemble tree-ring chronologies to reconstruct Northern Hemisphere summer temperatures over the Common Era.  The paper does not provide the chronologies and nor were data archived. The paper does (in the supplemental materials) describe how the 16 member site chronologies were developed.  This R script uses dplR and the information in Buntgen et al. (2020) with the raw ring width files provided by Buntgen et al. (2021) to create the 9 chronologies used in the reconstruction. 

The Tucson raw ring width files were archived here: https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/buentgen2021/

Because the chronologies used by Buntgen et al. (2020) are not available, the output of this script has not yet been able to be checked against those used in the paper. 

### References

Büntgen, U., Arseneault, D., Boucher, É., Churakova, O.V., Gennaretti, F., Crivellaro, A., Hughes, M.K., Kirdyanov, A.V., Klippel, L., Krusic, P.J. and Linderholm, H.W., 2020. Prominent role of volcanism in Common Era climate variability and human history. *Dendrochronologia*, 64, 125757.

Büntgen, U., Allen, K., Anchukaitis, K.J., Arseneault, D., Boucher, É., Bräuning, A., Chatterjee, S., Cherubini, P., Churakova, O.V., Corona, C. and Gennaretti, F., 2021. The influence of decision-making in tree ring-based climate reconstructions. *Nature Communications*, 12(1), 3411.
