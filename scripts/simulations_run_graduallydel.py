# vim: fdm=indent
'''
author:     Fabio Zanini
date:       16/08/13
content:    Run the graduallydel simulaions for fig 4A.
'''
# Modules



# Conditions
high = 0.003
ok = 0.001
low = 0.0003

script_fixation_loss_simulations_final.py ../data/simulations/gradually_del/ 10 1e4 1e-2 0.001 0.004 1.0 0.01 1 1 2e-5 0.01 0 0.1 0.2 0.2 0.4 0.4 0.6 0.6 0.9
