# Network method of cuts

## Goal of this project

This code enables the application of the network method of cuts as described in the paper: "Tilg, Ambühl, Batista, Menendez, Leclercq, Busch. From Corridor to Network Macroscopic Fundamental Diagrams: A Semi-analytical Approximation Approach." It is implemented for the case study of the Sioux Falls network. The repo includes the data from a Cell Transmission Model ground truth as well as from the applied state-of-the-art methods.

## Running the code
The code was developed and tested with MATLAB 2023a.

To apply the method, the `src/main.m` file needs to be run. The input data needed is mainly the network. It needs to be specified as the `src/data/network_1.m` file. The fundamental diagram is specified in the `src/main.m`. The specification of the network variational theory can be found in `src/functions/hypernetwork/nvt.m`.

## References

The mathematical background of the methods implemented in this repository are described in:

- Tilg, Ambühl, Batista, Menenden, Leclercq, Busch. From Corridor to Network Macroscopic Fundamental
Diagrams: A Semi-analytical Approximation Approach.

The "network variational theory" methods, i.e. the nvt.m function, was developed in:

- Tilg, G., Ambühl, L., Batista, S., Menendez, M., & Busch, F. (2021). On the application of variational theory to urban networks. Transportation Research Part B: Methodological, 150, 435-456.

The function MFD.m is borrowed from the code for the original method of cuts as described in:

- Leclercq, L., & Geroliminis, N. (2013). Estimating MFDs in simple networks with route choice. Transportation Research Part B: Methodological, 57, 468-484.

If you use the respective codes, please refer to these papers.
