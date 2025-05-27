# Spike-contrast: A Measure of Spike Train Synchrony

Spike-contrast is a novel, time-scale-independent, and multivariate measure of spike train synchrony, as described in the following peer-reviewed publication:

**Ciba, M., Isomura, T., Jimbo, Y., Bahmer, A., & Thielemann, C. (2018).**
*Spike-contrast: A novel time-scale independent and multivariate measure of spike train synchrony.*
Journal of Neuroscience Methods, 293, 136–143.
DOI: [10.1016/j.jneumeth.2017.09.008](https://doi.org/10.1016/j.jneumeth.2017.09.008)
[Link to the publication](https://www.sciencedirect.com/science/article/abs/pii/S0165027017303370)

## About Spike-contrast

Spike-contrast was developed to address challenges in measuring spike train synchrony, offering a robust tool that is independent of time scale and capable of analyzing multivariate spike data. It provides a single synchrony value but also a synchrony curve (synchrony as function of the time scale). The measure has applications in neuroscience, including studies of neural coding and network dynamics. 

## Repository Overview

This repository hosts the MATLAB and Python implementations of the Spike-contrast algorithm. Below, we provide details about the respective implementations.

### MATLAB Implementation

* **Description**: This is the original MATLAB implementation used in the publication.
* **Status**: This code is no longer actively maintained or updated.
* **Usage**: Researchers may use this version for reproducing the results presented in the publication.
* **Limitations**: For new applications, consider using the up-to-date Python implementation described below.

### Python Implementation

* **Description**: This is the original Python implementation of Spike-contrast.
* **Status**: While this version is functional, it may be outdated.
* **Current Version**: The most up-to-date Python implementation is integrated into the `Elephant` Python library.
* **Documentation**: Detailed documentation and usage examples are available in Elephant's official documentation.

  * [Spike-contrast in Elephant](https://elephant.readthedocs.io/en/latest/reference/_toctree/spike_train_synchrony/elephant.spike_train_synchrony.spike_contrast.html)

## Getting Started

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/yourusername/spike-contrast.git
   cd spike-contrast
   ```

2. **MATLAB Usage**:
   Open the `.m` files in MATLAB and follow the provided examples.

3. **Python Usage**:
   Run the `.py` files in your Python environment, or use the Elephant library for the latest implementation.

   To install Elephant:

   ```bash
   pip install elephant
   ```

4. **Further Information**:
   Refer to the publication for theoretical details and performance analysis of the Spike-contrast measure.

## License

This project is distributed under the MIT License. See the `LICENSE` file for more details.

## Citation

If you use Spike-contrast in your research, please cite the original publication:
*Ciba, M., Isomura, T., Jimbo, Y., Bahmer, A., & Thielemann, C. (2018). Spike-contrast: A novel time scale independent and multivariate measure of spike train synchrony. Journal of neuroscience methods, 293, 136-143.*
[DOI: 10.1016/j.jneumeth.2017.09.008](https://doi.org/10.1016/j.jneumeth.2017.09.008)

## Contact

For inquiries or issues related to this repository, please contact:

Manuel Ciba

manuel.ciba@th-ab.de

or

BioMEMS lab
(University of Applied Sciences Aschaffenburg)

biomems.lab.thab@gmail.com


