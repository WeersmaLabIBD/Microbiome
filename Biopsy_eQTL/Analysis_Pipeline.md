---
title: "eQTL analysis based on mucosal biopsy RNA-seq from patients with IBD"
creator: "Shixian"
date: "03/14/2019"
---

# eQTL analysis based on mucosal biopsy RNA-seq from patients with IBD

This project is to identify the eQTL effect in context of inflammation and non-inflammation in mucosal biopsy in IBD
-----------------------------
---
RNA-seq Data: "171 individuals; 299 biopsy"
Genomic Data: "171 individuals; WES+GSA"
Sample Excluded: "8"
Sample Included: "185CD + 106UC+IBDU"
---

Models used:
 - Model 1 (simple fixed model)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝜀
```
 - Model 2 (add random effect)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠 + 𝜀
```
 - Model 3 (add interaction term between SNPs and inflammation)
```
𝑔𝑒𝑛𝑒 𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛=𝛼 + 𝛽𝑆𝑁𝑃 + 20𝛽𝑒𝑥𝑝𝑟𝑒𝑠𝑠𝑖𝑜𝑛𝑃𝐶𝑠 + 𝑟𝑒𝑙𝑎𝑡𝑒𝑑𝑛𝑒𝑠𝑠 + 𝛽𝑆𝑁𝑃×𝑖𝑛𝑓𝑙𝑎𝑚𝑚𝑎𝑡𝑖𝑜𝑛 + 𝜀
```
