# BacktestingEV


***If you need to back-test your forecasting model as complexe as it might be, you have stumbled upon the solution !***


## Installation

```
pkg> add https://github.com/AMauro130/BacktestingEV.jl.git
```
Then
```
using BacktestingEV
```


## Presentation

This package will allow you to **back-test your forecasting model**.
We provide you a certain amount of functions to be applied to a prediction model within the known data.

You will be able to use some simple back-testing functions (as analysing the average values of prediction) and also more complexe ones to evaluate the extreme values of the model.
A meaninful part of the work has been granted to the **Extreme Values Theory (EVT)**.

*Remark : throughout this project we note Extreme Values as EV.*


## Requirement

To use this package properly, you need to provide _2 objects_ : **1 vector** (1 by n) for the *real data* and **1 matrix** (m by n) containing the predictions.

*Note : (n rows by m columns)*


## How to choose the appropriate threshold selection algorithm

**This table helps you to choose the algorithm**

To extract the extreme values from a distribution, you can choose from 6 algorithms based on very different methods.


<div id="header" align="center">
	<img src="https://user-images.githubusercontent.com/92920225/180974919-b05b1df7-ec06-45cf-812f-794a0ccb2595.png" alt="help table">
</div>

