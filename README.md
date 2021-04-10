# Single Electron Transistor

## Introduction
This repository shows the characteristics of a Single Electron Transistor (SET). In this device the electrons flow through a tunnel junction between source/drain to a quantum dot (conductive island). Moreover, the electrical potential of the island can be tuned by a third electrode, known as the gate, which is capacitively coupled to the island. The conductive island is sandwiched between two tunnel junctions, which are modeled by a capacitor C<sub>d</sub> and C<sub>s</sub> and a resistor R<sub>d</sub> and R<sub>s</sub> in parallel. 

## The code
* The functions [f_tunnel.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/f_tunnel.m) and [f_tunnel0.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/f_tunnel0.m) compute the tunneling rate

* The functions [master_equation.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/master_equation.m) and [master_equation0.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/master_equation.m) compute the master equation

* The function [monte_carlo.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/monte_carlo.m) calculates the Monte Carlo method

For more information type:
* `help f_tunnel`
* `help f_tunnel0`
* `help master_equation`  
* `help master_equation0`
* `help monte_carlo` in Matlab.

## How to run
Open Matlab and run [SET_characteristics.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/SET_characteristics.m) to see some characteristics of a single electron transitor or run [montecarlo_masterequation_comparison.m](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/montecarlo_masterequation_comparison.m) to see the comparison between Monte Carlo and Master Equation method to calculate current in a SET.

There is also a GUI [SET_characteristics_app.mlapp](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/SET_characteristics_app.mlapp) made up to easy understand the SET characteristics.

Application design

![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/Set_Behaviour/SET_app.png)

## Results

### SET characteristics

Drain current vs. drain voltage
![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/Set_Behaviour/drain_current_vs_drain_voltage.png)

Electrons in dot
![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/Set_Behaviour/electrons_in_dot.png)

dI/dV vs. drain voltage
![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/Set_Behaviour/dI_dV_vs_drain_voltage.png)

Coulomb blockade
![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/Set_Behaviour/coulomb_blockade.png)

### Comparison between Monte Carlo and master equation

Tunneling rate
![](https://github.com/MatteoOrlandini/SET-Simulation/blob/master/Images/MonteCarlo_MasterEquation/Tunneling_rate.png)

Charge density for different drain voltages
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/charge_density_for_different_drain_voltages.png)

Drain current vs. drain voltage at fixed gate voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/current_vs_drain_voltage_fixed_gate_voltage.png)

Charge density for different gate voltages
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/charge_density_for_different_gate_voltages.png)

Drain current vs. gate voltage at fixed drain voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/current_vs_gate_voltage_fixed_drain_voltage.png)

Sum of the number of charges that arrived and left the drain
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/sum_of_the_number_of_charges_arrived_and_left_drain.png)

Number of charges arrived at the drain at fixed gate voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/number_of_charge_arrived_at_drain_fixed_gate_voltage.png)

Number of charges that left at the drain at fixed gate voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/number_of_charge_that_left_drain_fixed_gate_voltage.png)

Drain current vs. time at fixed drain voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/drain_current_vs_time_fixed_drain_voltage.png)

Drain current vs. time at fixed gate voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/drain_current_vs_time_fixed_gate_voltage.png)

Drain current vs. drain voltage at fixed gate voltage comparison between Monte Carlo and master equation
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/drain_current_vs_drain_voltage_fixed_gate_voltage_comparison.png)

Drain current vs. gate voltage at fixed drain voltage comparison between Monte Carlo and master equation
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/drain_current_vs_gate_voltage_fixed_drain_voltage_comparison.png)

Charges stored into the dot varying gate and drain voltage
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/charges_stored_into_dot_varying_gate_drain_voltage.png)

Drain current varying gate and drain voltage comparison between Monte Carlo and master equation
![](https://github.com/MatteoOrlandini/Single-Electron-Transistor/blob/master/Images/MonteCarlo_MasterEquation/drain_current_varying_gate_drain_voltage_comparison.png)
