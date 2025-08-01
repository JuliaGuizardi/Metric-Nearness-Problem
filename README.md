The Metric Nearness Problem consists of finding the closest distance matrix, using the $\ell_p$ norm, to a given dissimilarity matrix such that metric properties, particularly the triangle inequalities, are satisfied. This problem can be applied in various contexts, including image processing, data clustering, and sensor location (when additional constraints are imposed), where noisy or incomplete distance data must be corrected.  
Two numerical optimization methods are studied and applied: the Augmented Lagrangian method via Algencan, and the Dykstra projection method, which was specifically programmed for this problem.

## Algencan

In the `Algencan` folder, the user will find the file `algencan_l2.f90`, which contains the code to solve the metric nearness problem with the $l_2$ norm using ALGENCAN.
To use this code, it is necessary to download ALGENCAN from the [TANGO website](https://www.ime.usp.br/~egbirgin/tango/). The version used for this project was **ALGENCAN 3.1.1**, and other versions may not work.

1. Compile the file by typing:

   ```
   gfortran -O3 algencan_l2.f90 -L$ALGENCAN/lib -lalgencan -lhsl -o algencan
   ```
   

2. Run by typing:

   ```
   ./algencan
   ```

3. A prompt will appear asking for the data file. Use `Examples/test_100.txt` as an example.


4. After the execution ends, another prompt will ask for the output file name.


The user can also define the input and output files directly in the code by setting:

```fortran
inputfnm = 'test_100.txt'
outputfnm = 'output_file.txt'
```

and commenting out the prompts.

## Dykstra

In the `Dykstra` folder, one will find code for solving the problem using the \$l\_1\$, \$l\_2\$, and \$l\_\infty\$ norms.

1. Compile by typing:

   ```bash
   gfortran -O3 dykstra_l2.f90 -o dykstra_l2
   ```

2. Run by typing:

   ```bash
   ./dykstra_l2
   ```

3. A prompt will appear asking for the data file. Use `Examples/test_100.txt` as an example.


4. After execution, a prompt will ask for the output file name.


As with Algencan, the user can define the input and output files directly in the code by setting:

```fortran
inputfnm = 'test_100.txt'
outputfnm = 'output_file.txt'
```

and commenting out the prompts.

To run Dykstra with the \$l\_1\$ or \$l\_\infty\$ norm, it is important to adjust **gamma**.
 A larger gamma leads to slower but more accurate convergence, while a smaller gamma results in faster (and sometimes less accurate) results. This value can be modified directly in the code.

## Examples

In the `Examples` folder, there are Julia language codes showing how the synthetic data was created and how real data was converted to be used in this problem.

The Fortran codes (Algencan and Dykstra) expect the matrices to be in a specific format. To transform any dissimilarity matrix into this format, you can use the ```Examples/create_data.ipynb``` notebook by calling the function:

```julia
function write_noisy_matrix_to_file(n, X)
```

where `n` is the size of the matrix and `X` is the dissimilarity matrix.

More about the project can be found in the Master’s Thesis: [metricnearnessproblem.pdf](https://github.com/JuliaGuizardi/Metric-Nearness-Problem/blob/main/metric-nearness-problem.pdf).

If you have any questions, feel free to contact me at **[ju\_guizardi@hotmail.com](mailto:ju_guizardi@hotmail.com)**.
