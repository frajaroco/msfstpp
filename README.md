# msfstpp: Mark summary functions for spatio-temporal point processes

Basic toolkit for the exploration and analysis of the spatio-temporal point patterns through mark summary functions. This repository is based on `spatstat`, `splancs`, `stpp`, `stats` and `KernSmooth` packages.

## Installation guide

The easiest way to install the development version of msfstpp from github is using the devtools package which can be installed run the next command:
```
install.packages('devtools', dependencies=TRUE)
```
and thereafter run the commands:
```
require(devtools)
devtools::install_github('frajaroco/msfstpp')
```

## References
- [Stoyan, D., Rodríguez-Cortés, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. *Spatial Statistics*. DOI: http://dx.doi.org/10.1016/j.spasta.2017.02.006.](http://www.sciencedirect.com/science/article/pii/S2211675317300696)
- [González, J. A., Rodríguez-Cortés, F. J., Cronie, O. and Mateu, J. (2016). Spatio-temporal point process statistics: a review. *Spatial Statiscts*, **18**:505-544.](http://www.sciencedirect.com/science/article/pii/S2211675316301130)
- Rodríguez-Cortés, F. J., Mateu, J. and Stoyan, D. (2017). Mark variograms and other mark-based second-order characteristics in the statistics of spatio-temporal point processes.

## CiteBibtex
```
@misc{r16,
	author = {Francisco J. Rodr\'{i}guez-Cort\'{e}s},
	title = {msfstpp: Mark summary functions for spatio-temporal point processes},
	year = {2016},
	note = {GitHub repository},
	url = {\url{https://github.com/frajaroco/msfstpp}}
}

```
