## ProSeqTeleporter: The Rapid Sequence Explorer

[![CI Build](https://github.com/bayer-int/proseq-teleporter/actions/workflows/build.yml/badge.svg)](https://github.com/bayer-int/proseq-teleporter/actions/workflows/build.yml)
[![Static Typing](https://github.com/bayer-int/proseq-teleporter/actions/workflows/static_typing.yml/badge.svg)](https://github.com/bayer-int/proseq-teleporter/actions/workflows/static_typing.yml)
[![Demo Page](https://github.com/bayer-int/proseq-teleporter/actions/workflows/generate_demo.yml/badge.svg)](https://github.com/bayer-int/proseq-teleporter/actions/workflows/generate_demo.yml)

### Setup environment
Compatible python version: 3.9, 3.10, 3.11

1. ```pip install --upgrade pip```
2. ```pip install -r requirements-dev.txt```

### Demo
You can view the demo notebook [here](https://legendary-adventure-oz4ywvn.pages.github.io/)
### Abstract
In the highly dynamic field of pharmaceutical R&D, the development of therapeutic biologics demands innovative solutions that enhance efficiency and reduce costs while maintaining robustness and reliability. We present ProSeqTeleporter, a state-of-the-art tool designed to address these needs by optimizing the process of multi-site mutagenesis, a crucial step in protein engineering.
ProSeqTeleporter empowers the creation of any desired mutant combinations from multiple positions of interest, with numerous variations for each position. It intelligently divides sequences with mutations, allowing for their reuse across multiple design-build-test-learn cycles, thereby significantly accelerating the construction of protein engineering libraries.
To illustrate the power of ProSeqTeleporter, consider the scenario of targeting any preferred combinations of mutants from 24 positions of interest, with 2 variations for each position. This allows for the creation of over 10^7 distinct combinations in the sequence space. The tool is designed to facilitate the one-step construction of any desired combination found within over 10^7 "teleportable" coordinates in the sequence space. This capability enables ProSeqTeleporter to move beyond the limitations of traditional protein library design. It provides a unique capability to instantly navigate, explore, and sample extensive sequence spaces in a single step, all while optimizing costs and ensuring a reliable process.
ProSeqTeleporter is more than just a multi-site mutagenesis tool. Its potential for seamless integration with machine learning processes further elevates its significance, positioning it as a key player in next-generation protein engineering strategies.
In summary, ProSeqTeleporter emerges as a valuable asset in the field of therapeutic biologics engineering. By enhancing the efficiency of the protein engineering process, it brings us one step closer to our goal: improving patient lives through the development of effective and affordable biological therapies.

![concept_picture.png](concept_picture.png)

### Process Overview
![ProcessOverview.png](ProcessOverview.png)


### Data Sources
- **Codon usage:**\
https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10029

- **Ligation fidelity:**\
Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly
Vladimir Potapov, Jennifer L. Ong, Rebecca B. Kucera, Bradley W. Langhorst, Katharina Bilotti, John M. Pryor, Eric J. Cantor, Barry Canton, Thomas F. Knight, Thomas C. Evans Jr., and Gregory J. S. Lohman
ACS Synthetic Biology 2018 7 (11), 2665-2674
DOI: 10.1021/acssynbio.8b00333