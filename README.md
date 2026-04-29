# Coastal_cleanup_prioritisation
Here we provide the code for the GIS based beach cleanup prioritisation study.

Code files are numbered in the order in which they should be run/processed. If implementing these methods into a new study region, ensure you read in the function script, and run through the following steps.

1. Ensure your clean-up data provides at minimum enough categorial resolution to seperate debris that is easily dispersed over water and those that are not, ideally your clean up locations will be revisited over several years.
2. Calculate gravity for clean up sites using an spatially explicit population density layer such as worldpop (https://www.worldpop.org/). If your coastal location also includes island beaches, ensure you calculate distnace over ocean to nearest ports (as shown in the code). We used distance to roads to account for accessibility in over land travel, this is an optional component, and can be done with other spatial layers (such as distance to urban centres).
3. Calculate exposure for clean up sites using windfetch package and customised functions, steps are in the code.
4. It is important to train a model that identifies the relationship between gravity (population pressure) and exposure with debris accumulation in your area before predicting across your coastline, as these relationships will differ locally depending on your study area. See code for steps
5. Create coastal points, or use previously generated coastal points and calculate both gravity and expsosure for each point using the same methods as steps 2 and 3.
6. Apply the training model you created in step 4 to your coastal points to generate predicted debris accumulation (as outlined in the code)
7. Prioiritisation combining amenity value (using gravity) and debris accumulation values can be done either in R through bivariatemap package, or it can be generated in a GIS (such as arc or qgis) using a bivariate function (quartiles)

Don't hesitate to send an email if you have any questions about implementing the methods to your system!
