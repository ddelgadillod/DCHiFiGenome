ALLHiC Scaffolding
=

We used ALLHiC pipeline for the scaffolding, we used several strategies by varing the input data and reference genome used for the synteny:

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

- Strategy 0: 
    Reference genome for sytenty -RH89-039-16 v3 potato assembly as reference, Input data: Canu draft assembly v1 and raw HiC reads.
- Strategy 1: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and raw HiC reads.
- Strategy 2: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and HiC valid contacts.
- Strategy 3: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v1 and HiC valid contacts.
- Strategy 4: Phasing and scaffolding with RH89-039-16 v3 potato assembly as reference, Canu draft assembly v3(corrected with HiC valid contacts) and HiC valid contacts.
- Strategy 5: Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v1 and HiC valid contacts.

### Strategy 6
Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v2(corrected with raw HiC reads) and HiC valid contacts.

### Strategy 7
Phasing and scaffolding with DMv6 potato assembly as reference, Canu draft assembly v3(corrected with HiC valid contacts) and HiC valid contacts


