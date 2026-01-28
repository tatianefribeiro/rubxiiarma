# RUBXII-ARMA Model

This repository provides R scripts to reproduce the application presented in the paper
“A dynamical regression model for double-bounded time series based on the reflected unit Burr XII distribution.”

The analysis focuses on the monthly proportion of stored hydroelectric energy in Northern Brazil, as described in Section 7 of the article.

## Usage

Run the main script:

application.r – data analysis, model fitting, diagnostics, and forecasting


## Competing Models

Implementations of the βARMA and KARMA models were obtained from:

https://github.com/vscher/barma

http://www.ufsm.br/bayer/boot-barma.zip

https://github.com/fabiobayer/KARMA

Details are provided in application.r.


## Reference

If you use this code, please cite:

@article{Ribeiro_PenaRamirez_Guerra_Alencar_Cordeiro_RUBXIIARMA,
  author  = {Tatiane Fontana Ribeiro and Fernando A. Peña-Ramírez and Renata Rojas Guerra and Airlane P. Alencar and Gauss Moutinho Cordeiro},
  title   = {A dynamical regression model for double-bounded time series based on the reflected unit Burr XII distribution},
  journal = {Environmental and Ecological Statistics},
  year    = {2026}
}
