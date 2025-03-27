Traffic eaction Model (TRM) Stochasztikus Szimulációs Eszköz

Ez a MATLAB szimuláció egy többcsomópontos közlekedési úthálózat dinamikájának modellezésére szolgál, 
amely lehetővé teszi mind a sztochasztikus (Gillespie-algoritmus), mind a determinisztikus (ODE) megközelítés vizsgálatát.

Fájlok:
TRM_red_ode_mod.m: Determinisztikus ODE modell az nyílt hurkú TRM modell szimulációjához külső be- és kifolyásokkal
TRM_external_flows.m: Függvény az időben változó külső be- és kifolyások definiálásához
gillespie_simulation.m: Sztochasztikus szimuláció Gillespie-algoritmussal
main.m: Fő szkript a szimuláció futtatásához