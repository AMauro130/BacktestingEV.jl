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

To use this package properly, you need to provide _2 objects_ : **1 vector** (with length n) for the *real data* and **1 matrix** (m by n) containing the predictions.

*Note : (n rows by m columns)*


###### Example

From a CSV file :

```
data = convert_to_matrix("D:\\users\\your_project\\simulation\\data.csv")
model = convert_to_matrix("D:\\users\\your_project\\simulation\\model.csv")
```

*Here it is n = 1800 and k = 10000*

*Data*
<div id="header" align="center">
	<img src="https://user-images.githubusercontent.com/92920225/181226632-66a8719f-2f97-49a9-a1e3-b048c56bf298.png" alt="data example" width=20% height=20%>
</div>

*Model*
<div id="header" align="center">
	<img src="https://user-images.githubusercontent.com/92920225/181226837-36d93d00-334e-4fef-bca1-77ca70c1f9b1.png" alt="model	 example">
</div>


## How to choose the appropriate threshold selection algorithm

**This table helps you to choose the algorithm**

To extract the extreme values from a distribution, you can choose from 6 algorithms based on very different methods.


<div id="header" align="center">
	<img src="https://user-images.githubusercontent.com/92920225/180974919-b05b1df7-ec06-45cf-812f-794a0ccb2595.png" alt="help table">
</div>

