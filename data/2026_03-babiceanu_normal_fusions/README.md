From [Recurrent chimeric fusion RNAs in non-cancer tissues and cells](https://doi.org/10.1093/nar/gkw032)

recurrent_normal.clean.csv was modified 

- replace non-delimiter "-" with "."
    - example: VAMP1-CD27-AS1 -> VAMP1-CD27.AS1
- removed extra spaces after fusions
- removed a duplicate fusion rows
    - keep one with greater recurrency reported
        - *SCNN1A-TNFRSF1A ,9,5
        - SCNN1A-TNFRSF1A ,7,5
        - *TIMM23B-LINC00843,51,26
        - TIMM23B-LINC00843 ,7,6
        - *DHRS1-RABGGTA,13,8
        - DHRS1-RABGGTA,9,5
