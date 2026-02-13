# Adolescence Self-Other Transfer

**Repository for: Kingston, Ellett, Richards, Burgess, & Barnby (2025).**

## Project Overview
This repository contains the data and computational code for a study investigating how humans generalize information from themselves to others (**self-to-other transfer**) and from others to themselves (**other-to-self transfer**). Focusing on the sensitive period of adolescence, this research explores the mechanisms of social generalization and how they underpin the development of stable social representations.

## Theoretical Framework:
The study is built upon two core mechanisms of social cognitive theory:
* **Self-Insertion (Social Projection):** The process where a "self" initially assumes others to be like their own representation, using their own preferences as a prior for predicting others.
* **Social Contagion:** The mechanism through which exposure to a partner’s preferences elicits a shift in one's own sense of self, leading toward an "interpersonal steady state".

## Experimental Paradigm: The Intentions Game
The task consists of three distinct phases[cite: 701]:
1. [**Phase 1 (Self-Baseline):** Participants make forced choices between options to establish their baseline social preferences ($\theta_{ppt}$).
2. **Phase 2 (Social Learning):** Participants observe and learn to predict the decisions of an anonymous partner, receiving feedback to update their partner model ($\theta_{par}$).
3. **Phase 3 (Post-Exposure):** Participants make choices for themselves again. This phase measures **Social Contagion** by assessing how the participant’s preferences shifted following the partner exposure.

## Repository Structure
* **`/Data`**: Contains behavioral data from the Intentions Game & Questionnaires. Data Dictionary Included.
* **`/Models`**: R and MATLAB scripts for Bayesian generative models[cite: 81].
* **`/HBI`**: Scripts for Hierarchical Bayesian Inference (HBI) used for concurrent model fitting and group-level comparison.
* **`/data_sim`**: Scripts for model-based simulations and agent-based modeling of group attractor dynamics.
* **`Beyond_ToM_Functions.R`**: Core utility functions for Bayesian updates, Fehr-Schmidt utility calculations, and belief convergence simulations.

## Citation
If you use this code or data, please cite:
> Kingston, J., Ellett, L., Richards, L., Burgess, H., & Barnby, J. M. (2026). *Social exclusion increases paranoia and reduces self and other learning flexibility in adolescents.* Clinical Psychological Science.

Related theoretical work:
> Barnby, J. M., Dayan, P., & Bell, V. (2023). *Formalising social representation to explain psychiatric symptoms*. Trends in Cognitive Sciences.

> Barnby, J. M., Nguyen, J., Griem, J., Wloszek, M., Burgess, H., Richards, L. J., ... & London Personality and Mood Disorders Consortium. (2025). *Self-other generalisation shapes social interaction and is disrupted in borderline personality disorder.* Elife, 14, RP104008.
---
**License:** This project is licensed under the AGPL-3.0 License.
