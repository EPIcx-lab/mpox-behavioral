# Code for "Adaptive behavior in response to the 2022 mpox epidemic in the Paris region" by D. Maniscalco et al.

## Description
This code serves to simulate infectious diseases, notably mpox, among men-who-have-sex-with-men (MSM) in the Paris region. The code simulates an agent-based process on the input temporal network

## Getting Started

### Dependencies
#### C++ Requirements
* Compiler: Apple Clang 14.0.3 (tested on macOS Ventura, ARM64)
* C++ standard: C++11 (required), but newer standards (C++14/17/20) are also supported
* No external libraries required (only standard C++ headers)

#### Python3 requirements
* Python 3.8.10
* Numpy 1.24.3
* Pandas 1.4.3

### Executing program

* Program is run with the command
```
python launcher.py
```
or 
```
python3 launcher.py
```
depending on the system. All the input parameters must be provided in the parameters.py file.

#### Input files

#### Output files
Each output file's name is a sequence of numbers separated by underscores. Numbers correspond to input parameters. The ordering of the numbers is guided in the 'launcher.py' file, and it is different

## Help
The code contains plenty of warning functions, that help to solve the most common problems and mistakes.

## Authors
Davide Maniscalco, davide.maniscalco@inserm.fr

## License
Maniscalco, D., Adaptive behavior in response to the 2022 mpox epidemic in the Paris region. Preprint at: https://www.medrxiv.org/content/10.1101/2024.10.25.24315987v1

Data and code belong to the authors:
Davide Maniscalco, Olivier Robineau, Pierre-Yves Boëlle, Mattia Mazzoli, Anne-Sophie Barret, Emilie Chazelle, Alexandra Mailles, Harold Noël, Arnaud Tarantola, Annie Velter, Laura Zanetti, Vittoria Colizza

If you use this material, please cite the above reference.

