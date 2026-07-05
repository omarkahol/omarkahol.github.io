---
layout: page
title: CMP++
description: A C++ library for Bayesian calibration of computer codes, with practical workflows for full, sequential, and modular inference.
img: assets/img/projects/cmp/logo.png
importance: 1
category: phd
---

<div class="repo p-2 text-center">
  <a href="https://github.com/omarkahol/cmp" rel="external nofollow noopener" target="_blank">
    <img class="repo-img-light w-100" alt="omarkahol/cmp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=cmp&amp;theme=default&amp;show_owner=false">
    <img class="repo-img-dark w-100" alt="omarkahol/cmp" src="https://github-readme-stats.vercel.app/api/pin/?username=omarkahol&amp;repo=cmp&amp;theme=dark&amp;show_owner=false">
  </a>
</div>

---

CMP++ is a C++ library used to perform the Bayesian calibration of deterministic computer models in the presence of measurement and model-form errors. It implements complete, sequential, and modular workflows so users can balance statistical rigor and computational cost depending on the application.

Approaches implemented include:
* Full Bayes
* Modular, sequential
* Modular, adaptive

## Documentation

For a deep dive into the mathematical framework, implementation details, and the API, please refer to the following resources:

* **Live Doxygen API Documentation**: [Live API Site](https://omarkahol.github.io/cmp/)

Below is a technical report that provides a comprehensive overview of the library's capabilities and usage:

<iframe src="https://docs.google.com/viewer?url=https://raw.githubusercontent.com/omarkahol/cmp/main/Technical_Doc/main.pdf&embedded=true" width="100%" height="800px" style="border: none; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);">
    <p>Your browser does not support viewing PDFs natively. <a href="https://github.com/omarkahol/cmp/blob/main/Technical_Doc/main.pdf">Click here to view or download the Technical Manual.</a></p>
</iframe>

## Getting Started

These instructions will help you compile and link the CMP++ library to your project. 

### Dependencies

The library depends on:
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [NLopt](https://nlopt.readthedocs.io/en/latest/)

### Installing

CMP++ can be compiled into a static library using the provided makefile. Make sure to modify the paths for the dependencies in your local environment.

    cd cmp
    make

## Running the Test Case

To compile and run the included tests, execute the following commands:

    cd tests
    make testname
    ./out_testname

*Note: In order to be able to view the generated plots and results, you must link against the matplotlibcpp library.*

## License

This project is licensed under the MIT License. See the [LICENSE.md](https://github.com/omarkahol/cmp/blob/main/LICENSE.md) file in the repository for details.

## Acknowledgments

This project has received funding from the European Union’s Horizon Europe research and innovation programme under the Marie Skłodowska-Curie grant agreement No 101072551 (TRACES). 

Please refer to the [TRACES project website](https://traces-project.eu) for additional details.