## To Do

- Get Argo data for SOCCOM float in region
- Check spectra to see if there is a scale dependence to the predictions
- Remove a simple mean to try and map only anomalies (for O2)
- Ask Lily about some details of 659 measurements.

## Discussion
- How to deal with error estimates for variables that are not measured by glider (like Nitrate):
  - train+test split. How to do this best? Maybe leave out profiles. If we leave out individual points then the surrounding data is very correlated.
  - Look at variance across decision trees
  - Leave one out cross validation
