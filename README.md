## Actin Dynamical Graph Grammar (DGG) Model of the Synaptic Spine Head
The code in this repository simulates membrane growth from an interaction between membrane and an internal actin cytoskeletal network for a synaptic spine head.

You can run and interact with the notebook itself using a version of Mathematica.

You will also need to install Plenum, the DGG modeling language package, in the right directory of your Wolfram installation.
It is in the directory Plenum_V1-22-26. After downloading it, follow the ReadMe.txt inside that folder.

Alternatively, you can download from https://computableplant.ics.uci.edu/theses/guy/Plenum.html
After download, in this current version of Plenum, change lines 606 and 610 of GrammarInterpreter.m, replacing both instances of "Random[]" with "RandomVariate[]" to allow multivariate thermal noise sampling.

As an example, a video of individual simulations is shown: circles being coarse-grained actin objects, color indicating angular energy with green being the barbed end, colored triangles being junction objects, squares being Arp2/3 branch objects, black triangles being capped-ends, and hexagons being CaMKII bundling objects.

For the data used to generate all plots in the manuscript (preprint: https://www.biorxiv.org/content/10.1101/2025.04.18.649558v2), download "SnapshotsInPaper:

Figure 2 is generated from snapshots in folders _STUB_AIP_1-9-26Arp, _STUB_AIP_1-9-26Cam, _STUB_AIP_1-9-26Cof, _STUB_AIP_1-9-26Aip1;

Figure 3 from snapshots in folders _STUB_AIP_2-4-26STUBMIXEDCameta1, _STUB_AIP_2-4-26STUBMIXEDCameta2, _STUB_AIP_2-4-26STUBMIXEDCameta4;

Figure 4 from snapshots in folders _STUB_AIP_1-11-26STUBMIXEDArpeta1, _STUB_AIP_1-11-26STUBMIXEDArpeta2, _STUB_AIP_1-11-26STUBMIXEDArpeta4;

Figure 5 from snapshots in folders _STUB_AIP_1-11-26STUBMIXEDCofeta1, _STUB_AIP_1-11-26STUBMIXEDCofeta2, _STUB_AIP_1-11-26STUBMIXEDCofeta4; and

Figure 6 from snapshots in folders _STUB_1-9-26Arp, _STUB_1-9-26Arp, _STUB_1-9-26Arp.

# 1-22-26 Version

![spinehead1](https://github.com/user-attachments/assets/4189c2b8-ee82-4855-a71d-15a308e350c4) ![spinehead2](https://github.com/user-attachments/assets/f6258df3-981a-4cf2-9b4f-2f8e44789154)


# 4-3-25 Version

![spinehead](https://github.com/user-attachments/assets/3ce5a407-11a6-4c15-af8e-160d8513172e)


2-19-26 Update

Additional work has been done for 3D version of the synaptic spine head, as outputted below:

<img width="356" height="403" alt="Screenshot 2026-02-19 at 10 56 13 PM" src="https://github.com/user-attachments/assets/d8e7b29f-07f2-4392-af5b-c02bc81d23f9" />

