# BacktestingEV



***If you need to back-test your forecasting model as complexe as it might be, you have stumbled upon the solution!***



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

You will be able to use some simple back-testing functions (as analysing the average values of prediction) and also more complexe ones to evaluate the extreme values of a model.
A meaninful part of the work has been granted to the **Extreme Values Theory (EVT)**.

*Remark : throughout this project we note Extreme Values as EV.*



## Requirements

To use this package properly, you need to provide _2 objects_ : **1 vector** (with length n) for the *real data* and **1 matrix** (m by n) containing the predictions.

(precise what kind of prediction array)

*Note : (n rows by m columns)*


###### Example

From a CSV file :

```
data = convert_to_matrix("D:\\users\\your_project\\simulation\\data.csv")
model = convert_to_matrix("D:\\users\\your_project\\simulation\\model.csv")
```

*Here it is n = 1800 and k = 10000*

*Data :*
<div align="left">
	<img src="https://user-images.githubusercontent.com/92920225/181226632-66a8719f-2f97-49a9-a1e3-b048c56bf298.png" alt="data example" width=20% height=20%>
</div>

*Model :*
<div align="center">
	<img src="https://user-images.githubusercontent.com/92920225/181226837-36d93d00-334e-4fef-bca1-77ca70c1f9b1.png" alt="model	 example">
</div>

You are now ready to back-test your model!


## How to choose the appropriate threshold selection algorithm

To extract the extreme values from a distribution, you can choose from **6 algorithms based on very different methods**.
From 8 datasets of distinct values with very particular distributions associated, we applied all these functions and compared the results with the expectations in the array below.

These expected values come from our visual intuition in terms of extreme values theory. This means, for example with this distribution :



Thus, the following table will help you to choose the most suitable algorithm (numbered from 1 ot 6)

<div align="center">
	<img src="https://user-images.githubusercontent.com/92920225/181231445-f95f5611-ec50-4c68-839a-86465a3a36f8.png" alt="help table">
</div>

<div align="center">
	These are percentage absolute relative error between (visual intuition and real return of the algorithms)
</div>


## References

Here are the links to the datasets used in the table for the threshold selecrtion :

(to fill in)