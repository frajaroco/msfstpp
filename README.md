# msfstpp: Mark summary functions for spatio-temporal point processes

Basic toolkit for the exploration and analysis of the spatio-temporal point patterns through mark summary functions. This repository is based on `splancs`, `spatstat`, `geoR`, `stpp`and `KernSmooth` packages.

## Installation guide

The easiest way to install the development version of `msfstpp` from GitHub is using the devtools package which can be installed run the next command:
```
install.packages('remotes')
```
and thereafter run the commands:
```
require(remotes)
install_github('frajaroco/msfstpp')
```
**warning:** If you have problems to load the repository. It is  is necessary to load the additional *R* package `plot3D` in order to plot the outputs of some functions.

## References
- Ballani, F., Rodríguez-Cortés, F. J., Mateu, J. and Stoyan, D. (2017). Mark-based second-order characteristics in the statistics for spatio-temporal point processes.
- [Stoyan, D., Rodríguez-Cortés, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. *Spatial Statistics*. **20**:125-147.](http://www.sciencedirect.com/science/article/pii/S2211675317300696)
- [González, J. A., Rodríguez-Cortés, F. J., Cronie, O. and Mateu, J. (2016). Spatio-temporal point process statistics: a review. *Spatial Statiscts*, **18**:505-544.](http://www.sciencedirect.com/science/article/pii/S2211675316301130)

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
