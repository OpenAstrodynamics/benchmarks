# ICATT 2016 Benchmarks

## Average Performance Relative to Fortran (N=100000)

| Problem      | Elements | Kepler | Lambert | Dopri |
| --           | --       | --     | --      | --    |
| C++          | 0.42     | 0.2    | 1.46    | 1.28  |
| Julia        | 2.29     | 0.71   | 1.1     | 2.97  |
| Java         | 7.54     | 0.7    | 3.2     | 19.22 |
| Python+Numba | 5.57     | N/A    | 1.56    | N/A   |
| Python       | 251.41   | 25.85  | 133.23  | N/A   |
| Matlab       | 186.41   | 282.96 | 196.89  | N/A   |

![Time vs. SLOC](avg_vs_sloc.svg)
