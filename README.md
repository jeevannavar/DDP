# MORONET: Multi-Omics gRaph cOnvolutional NETworks
Tongxin Wang, Wei Shao, Zhi Huang, Haixu Tang, Jie Zhang, Zhengming Ding, and Kun Huang

MORONET is a novel multi-omics data integrative analysis framework for classification tasks in biomedical applications.

![MORONET](https://github.com/txWang/MORONET/blob/master/MORONET.png "MORONET")
Overview of MORONET. \
<sup>MORONET combines GCN for multi-omics specific learning and VCDN for multi-omics integration. For concise illustration, an example of one patient is chosen to demonstrate the VCDN component for multi-omics integration after omics-specific learning. Pre-processing is first performed on each omics data type to remove noise and redundant features. Omics-specific GCN learns class prediction using omics features and the corresponding patient similarity network generated from the omics data. Cross-omics discovery tensor is calculated from initial predictions from GCN and forwarded to VCDN for final prediction. MORONET is an end-to-end model and all networks are trained jointly.<sup>

## Dependencies
* Python 3.6.5
* numpy 1.16.2
* pandas 0.23.0
* scikit-learn 0.21.1
* pyTorch 1.3.1

## Files
*main_brca.py*: An example of BRCA dataset for PAM50 subtype classification\
*models.py*: MORONET model\
*train_test.py*: Training and testing functions\
*load_data.py*: Read data from csv and perform data pre-processing\
*utils.py*: Supporting functions
