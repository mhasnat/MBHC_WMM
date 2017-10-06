# MBHC_WMM
Implementation of the clustering method: Model Based Hierarchical Clustering using Watson Mixture Model (MBHC-WMM)

- The MBHC-WMM method is an automatic method to cluster 3 dimensional axial data. This repo provides GUI demo with MATLAB code to do the following tasks: <br>
**a.** load 3 dimensional axial (labeled/unlabeled) data and display them. <br>
**b.** Generate/Synthesize 3 dimensional axial (labeled/unlabeled) data with given number of groups/cluster. <br>
**c.** Automatically cluster 3 dimensional axial data.  <br>

## How to use demo:
**Run** the MATLAB file name: **mbc\_wmm.m** <br>
**Load** data/samples files name: **S\_10000\_5\_Cl\_1.mat** or **S\_10000\_5\_Cl\_45.mat** <br>
_OR_ **Generate** 3D samples with specific number of cluster (edit text box - _Num Class_) <br>

## Application:
This clustering method has been used to cluster image normals (3D directional unit vectors) for analyzing depth/3D images. For details and other possible applications please see the reference. 

## Advantage:
This clustering method does not require to specify the number of clusters. It automatically determines it (see the reference).

## Limitations:
This clustering method is limited to cluster 3 Dimensional directional samples only.

## Extensions and scopes:
You can extend this method for higher dimensional samples.

## Reference:
1. Abul Hasnat, Olivier Alata and Alain Tremeau, Unsupervised Clustering of Depth Images using Watson Mixture Model, In Int. Conference on Pattern Recognition (**ICPR**) **2014**, Stockholm, Sweden.

2. Abul Hasnat, Olivier Alata and Alain Trmeau, Joint Color-Spatial-Directional clustering and Region Merging (JCSD-RM) for unsupervised RGB-D image segmentation, In IEEE Trans. on Pattern Analysis and Machine Intelligence (**TPAMI**), Vol 38, Issue 11, pages 2255 - 2268, **2015**.

