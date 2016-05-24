# Multi Expression Programming

## Source code

There are 6 source codes in this repository:

basic_mep.cpp - basic MEP algorithm. Solves regression and binary classification problems. This should be the starting point for everyone.

mep_multiple_pop.cpp - MEP with multiple (sub)populations for regression and binary classification problems.

mep_circuits.cpp - MEP for evolving digital circuits.

mep_threads.cpp - MEP with multiple (sub)populations evolved in different threads for regression and binary classification problems.

mep_multi_class.cpp - MEP algorithm. Solves regression and multi-class classification problems.

mep_parity_adfs.cpp - MEP algorithm with Automatically Defined Functions for even parity problems.

## How to run

Take a .cpp program from src folder and compile it.
You also need a dataset from the dataset folder. Make sure that you specify the path correctly.

## Datasets

There are 5 datasets in this repository (check data folder):

building1.txt - this is a symbolic regression problem taken from PROBEN1. It is used by mep.cpp, mep_multiple_pop.cpp, mep_multi_class.cpp and mep_threads.cpp for training.

cancer1.txt - this is a classification problem taken from PROBEN1. It is used by mep.cpp, mep_multiple_pop.cpp, mep_multi_class.cpp and mep_threads.cpp for training.

2x2_multiplier.txt, 3x3_multiplier.txt - this is a problem of designing digital circuits (for 2x2 and 3x3 multiplication). It is used by mep_circuits.cpp for training.

iris.txt - a classification problem with 3 classes. It is used by mep_multi_class.cpp for training.

gene1.dt - a classification problem with 3 classes. It is used by mep_multi_class.cpp for training.

even_6_parity.txt - even parity problems with 7 inputs. It is used by mep_parity_adfs.cpp for training.

#More info:

http://www.mep.cs.ubbcluj.ro

#Discussion Group:

https://groups.google.com/d/forum/mepx

#Contact author:

mihai.oltean@gmail.com
