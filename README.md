# eccomp2020
EC-Comp2020 problem: Designing a random number sequence to entertain game players. This problem is an N-variable, M-objective, L-constraint problem of creating cognitively unbiased random number sequences.

## Install
```
$ git clone https://github.com/opthub-org/eccomp2020.git
$ cd eccomp2020
$ pip install -r requirements.txt
```

## Usage
```
$ ./rngbias.py < rand50.txt
{"objective": 81.14285714285714, "constraint": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "error": null}
```

See help for detailed usage.
```
$ ./rngbias.py --help
```

For convenience, a random number generator is provided.
```
$ ./rand.py | ./rngbias.py
{"objective": 81.14285714285714, "constraint": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], "error": null}
```

See help for detailed usage.
```
$ ./rand.py --help
```

## Environmental Variables
Using environmental variables, you can specify the number of variables (i.e., the length of random number sequence), the definition of objective functions, and the parameters of cognitive bias.

|Variable              |Default                                                |
|----------------------|-------------------------------------------------------|
|`RNGBIAS_VARIABLES`   |`50`                                                   |
|`RNGBIAS_OBJECTIVES`  |`[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]`|
|`RNGBIAS_CONSTRAINTS` |`[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]`              |
|`RNGBIAS_LOWER_BOUNDS`|`[.1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1, .1]`     |
|`RNGBIAS_UPPER_BOUNDS`|`[.9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9, .9]`     |
|`RNGBIAS_BIAS_ALPHA`  |`[2, 2, 2, 2, 2, 27, 5, 0,  0, 1, 0, 0, 1, 0, 0]`      |
|`RNGBIAS_BIAS_BETA`   |`[5, 5, 5, 5, 5, 30, 8, 1,  0, 3, 0, 1, 2, 0, 0]`      |
|`RNGBIAS_BIAS_GAMMA`  |`[3, 3, 3, 3, 3,  1, 1, 3, 10, 4, 4, 4, 4, 4, 4]`      |

The objective functions can be configured through `RNGBIAS_OBJECTIVES`, which is a 2D array of indices of cognitive features F1--F15 where each 1D array in the 2D array represents an objective function as a sum of specified features. The default value defines a single objective function that sums all the 15 features, which is the setting of EC-Comp 2020 Single-Objective Track. You may change it to the setting of EC-Comp 2020 Multi-Objective Track as follows:
- for M=2, `RNGBIAS_OBJECTIVES=[[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15]]`;
- for M=3, `RNGBIAS_OBJECTIVES=[[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]`;
- for M=5, `RNGBIAS_OBJECTIVES=[[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]]`;
- for M=7, `RNGBIAS_OBJECTIVES=[[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14, 15]]`.

The parameters of cognitive bias can be specified via `RNGBIAS_ALPHA`, `RNGBIAS_BETA`, and `RNGBIAS_GAMMA` environmental variables, each of which should be a numerical array of length 15, corresponding to the lower bound, the upper bound, and the penalty weight of features F1--F15. The default values are adopted from ([Temsiririrkkul et al. 2014, Table II](https://dspace.jaist.ac.jp/dspace/bitstream/10119/12995/1/21068.pdf)).

## Configuration File
You can also specify environmental variables by a file. See `config.yml` for details.