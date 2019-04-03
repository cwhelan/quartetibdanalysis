# QuartetIBDAnalysis
Tools for analyzing IBD state and copy number variation in quartet sequencing experiments

QuartetIBDAnalysis contains a set of tools for computing the IBD sharing state between siblings in quartets based on genome-wide SNP calls, and for analyzing copy number variation in those quartets through the use of IBD state.

# Installation

Prequisites: Java 8

To clone the repo and build with gradle:

```
$ git clone https://github.com/cwhelan/quartetibdanalysis.git
$ cd quartetibdanalysis/
$ ./gradlew installDist
```

# Usage

Get a list of commands to run:

```
$ ./quartetibdanalysis --list
Using wrapper script /humgen/gsa-hpprojects/dev/cwhelan/quartetibdanalysis/build/install/QuartetIBDAnalysis/bin/QuartetIBDAnalysis
Running:
    /humgen/gsa-hpprojects/dev/cwhelan/quartetibdanalysis/build/install/QuartetIBDAnalysis/bin/QuartetIBDAnalysis --help
USAGE: QuartetIBDAnalysis <program name> [-h]

Available Programs:
--------------------------------------------------------------------------------------
Quartet IBD Analysis Program Group:              Tools to analyze IBD state and copy number variation in quartets
    AnnotateCopyNumberVariantsWIthIBDState       Annotate a copy number VCF with IBD state
    SiblingIBD                                   Compute IBD blocks between siblings
    TargetDispersedDuplication                   Compute DD location based on IBD

--------------------------------------------------------------------------------------
```
