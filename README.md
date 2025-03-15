
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsdmm

<!-- badges: start -->

[![gsdmm status
badge](https://till-tietz.r-universe.dev/badges/gsdmm)](https://till-tietz.r-universe.dev/gsdmm)
<!-- badges: end -->

------------------------------------------------------------------------

`gsdmm` implements short text classification via Dirichlet Mixture
Models proposed by [Yin and Wang
2014](https://www.semanticscholar.org/paper/A-dirichlet-multinomial-mixture-model-based-for-Yin-Wang/d03ca28403da15e75bc3e90c21eab44031257e80?p2df).
It provides a fast `c++` implementation and R interface for the Gibbs
sampler described in the paper. Specifically, `gsdmm` implements the
Likelihood function allowing for multiple occurrences of the same word
in a given text (EQ4).

**Benefits:**  

- very space and time efficient
- unlike LDA it requires only an upper bound on the number of clusters

**Development:**  

- I am planning to add a tuning function for the alpha and beta
  parameters of the gibbs sampler

## Installation

You can install the development version of gsdmm from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("till-tietz/gsdmm")
```

## Usage

Here is a minimal working example.

``` r
# we lemmatize and tokenize creating a list of character vector representing each text
text <- c(
  "Rockets are amazing.",
  "Witnessing a rocket in flight is a marvel of engineering.",
  "We should take a rocket to Mars.",
  "Rocket",
  "Have you ever seen a cat?",
  "Cats are fun.",
  "Your cat seems sweet.",
  "Cat"
) |>
  tolower() |>
  gsub(pattern = "[[:punct:] ]+", replacement = " ") |>
  textstem::lemmatize_strings() |>
  text2vec::word_tokenizer() |>
  lapply(function(i) i[!i %in% stopwords::stopwords()])

set.seed(42)

gsdmm::gsdmm(texts = text, n_iter = 100, n_clust = 5, alpha = 0.1, beta = 0.01, progress = FALSE)
#> $cluster
#> [1] 2 0 0 0 4 4 4 4
#> 
#> $distribution
#> 14 x 3 sparse Matrix of class "dgCMatrix"
#>          0 2 4
#> rocket   3 1 .
#> amaze    . 1 .
#> witness  1 . .
#> flight   1 . .
#> marvel   1 . .
#> engineer 1 . .
#> take     1 . .
#> mar      1 . .
#> ever     . . 1
#> see      . . 1
#> cat      . . 4
#> fun      . . 1
#> seem     . . 1
#> sweet    . . 1
```
