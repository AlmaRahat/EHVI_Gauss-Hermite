# Efficient Approximation of Expected Hypervolume Improvement using Gauss-Hermite Quadrature

This is a repository for the code to generate and use Gauss-Hermite approximation for EHVI calculations. 

The required Python libraries are: 
- [Numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [Pandas](https://pandas.pydata.org/)
- [evolalgos](https://pypi.org/project/evoalgos/)

To install the exact EHVI code, you need to go to "EHVI/EHVI_2D" and "EHVI/EHVI_3D" directories in turn, and issue the following command. 
```
$ g++ -Ofast -o EHVI -std=c++0x *.cc
```

An illustration on how to use the code is provided in the "notebooks/Illustration of running experiments.ipynb". You would need [Jupyter notebook](https://jupyter.org/) to run that. 

The data from the experiments ran for the paper is avaialble in the "data/experiments/" directory. For the meaning of the headers of that file, please refere to the illustrative notebook, in particular, the last paragraph where we introduced the keys to the labels. 
