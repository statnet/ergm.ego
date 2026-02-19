# Bibliography Management with Rdpack

This package now uses the Rdpack package to manage bibliographic references in the documentation. This provides a centralized way to maintain references across all documentation files.

## Structure

- **inst/REFERENCES.bib**: Central bibliography file containing all references in BibTeX format
- **R source files**: Use `\insertRef{key}{ergm.ego}` in roxygen2 comments to cite references
- **man pages**: Auto-generated from R source files using roxygen2

## How to Update References

### Adding a New Reference

1. Add the BibTeX entry to `inst/REFERENCES.bib`
2. Use `\insertRef{BibKey}{ergm.ego}` in the roxygen2 `@references` section of your R source file
3. Run `roxygen2::roxygenise()` to regenerate the man pages

### Updating an Existing Reference

1. Edit the entry in `inst/REFERENCES.bib`
2. Run `roxygen2::roxygenise()` to regenerate the man pages

## Example

In R source file:
```r
#' @references
#' 
#' \insertRef{KrMo2017}{ergm.ego}
```

In inst/REFERENCES.bib:
```bibtex
@Article{KrMo2017,
  author = {Pavel N. Krivitsky and Martina Morris},
  title = {Inference for social network models...},
  journal = {Annals of Applied Statistics},
  year = {2017},
  volume = {11},
  number = {1},
  pages = {427--455},
  doi = {10.1214/16-AOAS1010}
}
```

## Key References Updated

The following key reference has been updated with correct publication details:
- **KrBoMo2020**: Previously showed "to appear", now correctly shows year 2021, volume 64, pages 69-87

## Requirements

- Rdpack package (added to DESCRIPTION Imports and RdMacros)
- roxygen2 for documentation generation
