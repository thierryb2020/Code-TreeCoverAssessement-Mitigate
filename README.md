# CGIAR/Mitigate+ Project - Agricultural Production Analysis

## Overview

This code, developed by Paul Bostyn, is part of the CGIAR/Mitigate+ project aimed at analyzing agricultural production losses and gains based on tree cover levels. The code is divided into three main parts:

1. **R code1.R:**
   - Produces total agricultural production losses and losses per hectare for each tree cover level defined by the variable 'objectifalgo1'.

2. **R code2.R:**
   - Calculates the net present value, integrating all costs and benefits across a wide range of tree cover levels for each pattern defined in the variable 'algoobj'.

3. **R_agroforestry.R:**
   - Generates results related to agroforestry, including areas suitable for agroforestry and the associated production values.

## Usage

To use this code, you have two options:

1. **Using Pre-formatted Data:**
   - You can directly use the formatted data available in the file `df2All.xlsx`. Ensure that you have this file in the same directory as the code files.

2. **Downloading Data (Open Access):**
   - Alternatively, you can download the required data from freely accessible online sources. Ensure the data is in a compatible format to reproduce the results.

## Instructions

Follow these steps to run the code:

1. Open each R script (`R code1.R`, `R code2.R`, `R_agroforestry.R`) in your preferred R environment.
2. Set the appropriate values for variables such as 'objectifalgo1' and 'algoobj' as needed.
3. Execute the scripts in order to obtain the desired results.

## Dependencies

Make sure you have the necessary R packages installed. You can install them using the following commands:

```R
# Example for installing required packages
install.packages("package_name")
