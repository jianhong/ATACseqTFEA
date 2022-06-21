# Change in 0.99.5
  - update the license from 'GPL (>=2)' to 'GPL-3'
  - add additional sections with guidelines for bug report and contribution in
  readme.md.
  - split the class and method into three files: R/AllClasses.R,
  R/AllGenerics.R, and methods-TFEAresults.R.
  - simplify the `show` outputs for `TFEAresults` object
  - use loop count no more than `MAX_LOOPS` to control the nested loops
  - return `ggplot` object for the `plotES` function
  - break the `TFEA` function into multiple functions
  - update muliple documentations.

# Change in 0.99.4
  - add schematic diagram.

# Change in 0.99.3
  - replace `import` by `importFrom`.
  
# Change in 0.99.2
  - remove the import of `utils`.
  
# Change in 0.99.1
  - export function `reduceByPercentage`
  - fix the bug in vignette
  - import utils in DESCRIPTION.

# Change in 0.99.0
  - prepare for submission
