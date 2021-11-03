ALLHiC Scaffolding
=

We used ALLHiC pipeline for the scaffolding, we used several strategies by varing the input data and reference genome used for the synteny:

=
Initial assembly
-Canu v1: the HiCanu assembly 
-Canu v2: the output of ALLHIC_corrector giving in input Canu v1 and the raw Hi-C data 
-Canu v3: the output of ALLHIC_corrector giving in input Canu v1 and the Hi-C valid contacts 

=

| Strategy  | Reference genome | Initial assembly | Hi-C reads |
| ------------- | ------------- | ------------- | ------------- |
| st0  | RH89-039-16 v3   | Canu draft assembly v1 | raw HiC reads |
| st1  | RH89-039-16 v3   | Canu draft assembly v2(corrected with raw Hi-C reads) | raw Hi-C reads|
| st2  | RH89-039-16 v3   | Canu draft assembly v2(corrected with raw Hi-C reads) | Hi-C valid contacts|
| st3  | RH89-039-16 v3   | Canu draft assembly v1  | Hi-C valid contacts|
| st4  | RH89-039-16 v3   | Canu draft assembly v3(corrected with Hi-C valid contacts)  | Hi-C valid contacts|
| st5  | DMv6    | Canu draft assembly v1  | Hi-C valid contacts|
| st6  | DMv6    | Canu draft assembly v2(corrected with raw Hi-C reads) | Hi-C valid contacts|
| st7  | DMv6    | Canu draft assembly v3(corrected with Hi-C valid contacts) | Hi-C valid contacts|

