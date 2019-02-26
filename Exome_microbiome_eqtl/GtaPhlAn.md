### GraphLAN

```
://bitbucket.org/nsegata/graphlan/downloads/graphlan_commit_6ca8735.zip

Export:

export PATH=“Absolute path”/graphlan_commit_6ca8735/:$PATH

Modules:

module load Python(2.7.9 version)
module load Biopython
module load numpy
module load matplotlib

core document：
guide_1.txt (bacteria document)
annot_0.txt 
annot_1.txt

core commands：
graphlan_annotate.py --annot annot_0.txt guide.txt guide_1.xml
graphlan_annotate.py --annot annot_1.txt guide_1.xml guide_2.xml
graphlan.py guide_2.xml step_2.png --dpi 300 --size 7.5
```
