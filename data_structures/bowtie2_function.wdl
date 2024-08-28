## Bowtie2 accepts several function parameters.
## The first term is a function type. Available function types are:
## - C - constant
## - L - linear
## - S - square-root
## - G - natural logarithm
## The constant and coefficient types may be negative and/or floating point numbers.
##
## An example input JSON entry might look like:
## ```
## {
##     "interval_seed_substrings": {
##         "function_type": "S",
##         "constant": 1,
##         "coefficient": 0.50,
##     },
## }
## ```

version 1.1

# See the function documentation in bowtie2
# (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#setting-function-options)
struct Bowtie2Function {
    String function_type
    Float constant
    Float coefficient
}
