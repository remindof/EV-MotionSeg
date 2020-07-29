# EV-MotionSeg
MATLAB code for paper "Event-Based Motion Segmentation by Motion Compensation"

This is non-official code for paper https://arxiv.org/abs/1904.01293 by Timo Stoffregen, Guillermo Gallego, Tom Drummond, Lindsay Kleeman, Davide Scaramuzza.
The algorithm has been simplified and the code does not aim to reproduce the original paper exactly but to study the idea in the paper.

## Something different
- No details of  Dirac delta function approximation in original paper, thus I manually set the gradient of delta function, see function findGradDelta in "updateMotionParam.m"
- Only linear warp has been considered.

## Results
The images show two waving hands that are moving in opposite horizontal directions.

The probabilities of event clusters during the iterations:

![image](https://github.com/remindof/EV-MotionSeg/blob/master/results/iter1-3_clusters.png)

IWE (Image of Warpped Events) after 3 iterations:

![image](https://github.com/remindof/EV-MotionSeg/blob/master/results/iter%3D3_IWE.png)

The scenario with only 1 cluster but the parameter of cluster number was set to 2:

![image](https://github.com/remindof/EV-MotionSeg/blob/master/results/cluster_no%3D2but_only_1.png)
