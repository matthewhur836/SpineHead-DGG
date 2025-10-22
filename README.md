## Actin Dynamical Graph Grammar (DGG) Model of the Synaptic Spine Head
The code in this repository simulates membrane growth from an interaction between membrane and an internal actin cytoskeletal network for a synaptic spine head.

Otherwise, you can run and interact with the notebook itself using a version of Mathematica.

You will also need to install Plenum, the DGG package, in the right directory of your Wolfram installation: https://computableplant.ics.uci.edu/theses/guy/Plenum.html

After download, in this current version of Plenum, change lines 606 and 610 of GrammarInterpreter.m, replacing both instances of "Random[]" with "RandomVariate[]" to allow multivariate Hessian sampling.

As an example, a video of individual simulations is shown: circles being coarse-grained actin objects, color indicating angular energy with green being the barbed end, colored triangles being junction objects, squares being Arp2/3 branch objects, black triangles being capped-ends, and hexagons being CaMKII bundling objects.


# 10-21-25 Version 3

![out](https://github.com/user-attachments/assets/db7c5e77-196f-42a6-94fc-fc2f5ea2f12c)


# 4-3-25 Version 1

![spinehead](https://github.com/user-attachments/assets/3ce5a407-11a6-4c15-af8e-160d8513172e)

Preprint is available at https://www.biorxiv.org/content/10.1101/2025.04.18.649558v1

