# RPPF
Ranking of Pareto solutions based on projection free-energy


# Requirements
Python >= 3.7


# Usage

## Parameters

**rho: int** (parameter in augmented weighted Tchebycheff)

**num_top: int** (number of proposals)

**wind: float** (window of weights)

**Tstar: float** (value of Tstar)


## Target dataset

The target dataset is prepared as data.csv.

The example is (two-dimensional property space) as follows.

|  name  |  objective1  |  objective 2 |
| ------ | ------------ | ------------ |
|  Sc2O3 |  3.75        |  0.672       |
|  GeO2  |  3.05        |  0.518       |
|  :     |  :           | :            |

In present, two- and three-dimensional property space can be targetted.

## Execution
```
python rppf.py 
```

# License
This project is licensed under the terms of the MIT license.
