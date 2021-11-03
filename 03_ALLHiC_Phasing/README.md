ALLHiC Scaffolding
=

We used ALLHiC pipeline for the scaffolding, we used several strategies by varing the input data and reference genome used for the synteny:

----------------
#### Initial assembly:

- Canu v1 assembly: the HiCanu assembly 
- Canu v2 assembly: the output of ALLHIC_corrector giving in input Canu v1 and the raw Hi-C data 
- Canu v3 assembly: the output of ALLHIC_corrector giving in input Canu v1 and the Hi-C valid contacts 

----------------
#### Reference genome:

- [RH89-039-16](https://www.nature.com/articles/s41588-020-0699-x)
- [DMv6](https://academic.oup.com/gigascience/article/9/9/giaa100/5910251)

----------------

| Strategy  | Reference genome | Initial assembly | Hi-C reads |
| ------------- | ------------- | ------------- | ------------- |
| st0  | RH89-039-16   | Canu v1 assembly | raw HiC reads |
| st1  | RH89-039-16   | Canu v2 assembly | raw Hi-C reads|
| st2  | RH89-039-16   | Canu v2 assembly | Hi-C valid contacts|
| st3  | RH89-039-16   | Canu v1 assembly  | Hi-C valid contacts|
| st4  | RH89-039-16   | Canu v3 assembly  | Hi-C valid contacts|
| st5  | DMv6    | Canu v1 assembly  | Hi-C valid contacts|
| st6  | DMv6    | Canu v2 assembly | Hi-C valid contacts|
| st7  | DMv6    | Canu v3 assembly | Hi-C valid contacts|

