# ValCont

R package for content validity analysis

\strong{Package:} ValCont\
\strong{Title:} R Package for content validity coefficient estimation\
\strong{Version:} 0.1.0\
\strong{Date:} 07-28-2025\
[\strong{Authors\@R}](mailto:Authors@R)\strong{:}
- \strong{Cesar Merino-Soto (aut, cre, ctb, rev)}
- Email: [sikayax\@yahoo.com.ar](mailto:sikayax@yahoo.com.ar)
- ORCID: https://orcid.org/0000-0001-8049-7069
- \strong{Jose Livia-Segovia (aut, ctb)}
- Email: [joselivia\@gmail.com](mailto:joselivia@gmail.com)
- ORCID: https://orcid.org/0000-0003-2226-3349
- \strong{Diego Livia-Ortiz (aut, ctb, cre, rev)}
- Email: [diegolivia\@hotmail.com](mailto:diegolivia@hotmail.com)
- ORCID: https://orcid.org/0000-0002-2107-3140

\strong{Description:}\
`ValCont` is a dedicated content validity package in R. `ValCont` implement the computation of several coefficients used in content validity studies, with data usually obtained from selected participants such as expert judges or experiential judges. The coefficients calculated by ValCont are:

-   CVC (Hernandez-Nieto, 2002)
-   CVI (Matuza, 1977; Lynn, 1986)
-   CVIR (Polit et al., 2007)
-   CVR (Lawshe, 1975)
-   Psv (proportion of substantive agreementy) & Csv (coefficient of substantive validity); (Anderson & Gerbing, 1991)
-   V (Aiken, 1980, 1985)
-   MER (Mean of Expert Ratings; Penfield, & Miller, 2004)

Given that these coefficients can be conceptualized and/or are handled as proportions, with limits at 0.0 and 1.0, `ValCont` also implements asymmetric confidence intervals for each coefficient, appropriate for small samples and non-normally distributed data (Penfield & Giacobbi, 2004; Wilson, 1927).

Some functions were added to estimate other relevant aspects of the content validity analysis:

-   Difference between two independent content validity coefficients: Confidence intervals for difference (Merino-Soto, 2018) and standardized difference (Merino-Soto, 2023)

-   Ratio of two independent content validity coefficients and confidence Intervals (MOVER-R approach; Zou, Donner, & Qiu, 2025)

-   Homogeneity of ratings (H coefficient; Aiken, 1980, 1985).

-   Basic functions to make graphs of results are also implemented.

\strong{Install:}

You can install the development version of `ValCont` from GitHub using:
```R
if(!"devtools" %in% row.names(installed.packages())){
  install.packages("devtools")
}
devtools::install_github("Diegolivia/ValCont")
```
\strong{References:}

-   Aiken, L. R. (1980). Content validity and reliability of single items or questionnaires. Educational and Psychological Measurement, 40, 955-959. <https://doi.org/10.1177/001316448004000419>
-   Aiken, L. R. (1985). Three coefficients for analyzing the reliability and validity of ratings. Educational and Psychological Measurement, 45, 131-142. <https://doi.org/10.1177/0013164485451012>
-   Anderson, J. C., & Gerbing, D. W. (1991). Predicting the performance of measures in a confirmatory factor analysis with a pretest assessment of their substantive validities. Journal of Applied Psychology, 76(5), 732–740. <https://doi.org/10.1037/0021-9010.76.5.732>
-   Hernandez-Nieto, R. A. (2002). Contributions to Statistical Analysis. Merida, Venezuela: Universidad de Los Andes.
-   Lawshe, C. H. (1975). A quantitative approach to content validity. Personnel Psychology, 28,563-575
-   Lynn, M.R. (1986). Determination and quantification of content validity. Nursing Research, 35, 382–385.
-   Martuza, V.R. (1977). Applying norm-referenced and criterion-referenced measurement in education. Boston: Allyn & Bacon
-   Merino-Soto, C. (2023). Aiken’s V Coefficient: Differences in Content Validity Judgments. MHSalud: Revista En Ciencias Del Movimiento Humano Y Salud, 20(1), 1-10. <https://doi.org/10.15359/mhs.20-1.3>
-   Merino-Soto, C. (2018). Confidence interval for difference between coefficients of content validity (Aiken's V): a SPSS syntax. Anales de Psicología, 34(3), 587-590. <https://dx.doi.org/10.6018/analesps.34.3.283481>
-   Penfield, R. D. & Giacobbi, P. R., Jr. (2004) Applying a score confidence interval to Aiken’s item content-relevance index. Measurement in Physical Education and Exercise Science, 8(4), 213-225. <https://doi.org/> 10.1207/s15327841mpee0804_3
-   Penfield, R. D., & Miller, J. M. (2004). Improving Content Validation Studies Using an Asymmetric Confidence Interval for the Mean of Expert Ratings. Applied Measurement in Education, 17(4), 359–370. <https://doi.org/10.1207/s15324818ame1704_2>
-   Polit, D.F., Beck, C.T. and Owen, S.V. (2007), Is the CVI an acceptable indicator of content validity? Appraisal and recommendations. Res. Nurs. Health, 30: 459-467. <https://doi.org/10.1002/nur.20199>
-   Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association, 22, 209-212. <https://doi.org/10.2307/2276774>
-   Zou, G., Donner, A. and Qiu, S. (2025). MOVER-R for Confidence Intervals of Ratios. In Wiley StatsRef: Statistics Reference Online (eds N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch, F. Ruggeri and J.L. Teugels). <https://doi.org/10.1002/9781118445112.stat08085>

\strong{Depends:} R (\>= 2.10)\
\strong{Imports:} boots, ggplo2, stats, utils\
\strong{License:} GPL-3\
\strong{Encoding:} UTF-8\
\strong{LazyData:} true\
\strong{Maintainer:} Diego Livia-Ortiz\
\strong{URL:} https://github.com/Diegolivia/ValCont/

