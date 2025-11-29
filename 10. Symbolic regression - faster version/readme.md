# Symbolic regression with MEP

Multi Expression Programming - basic source code for solving symbolic regression problems

Faster version, which does not copy the values of variables in the *eval_matrix*.

*training_data* matrix is stored as an array of columns (columns are variables) instead of (naturally) an array of rows (rows are data)

The differences between this version and the one stored in the (01. Symbolic regression)[https://github.com/mepx/mep-basic-src/tree/master/01.%20Symbolic%20regression] folder are visible in the *compute_eval_matrix* and *fitness_regression* functions.


## Papers to read

Oltean Mihai, Dumitrescu D., [Multi Expression Programming](../papers/oltean_mep.pdf), Technical report, Babeș-Bolyai University, 2002.