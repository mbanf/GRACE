"""GRACE: gene regulatory network inference accuracy enhancement with Markov random fields.

Python implementation of Banf & Rhee, Scientific Reports 7:41174 (2017).

    from grace.pipeline import prepare, ensemble_network
    from grace.reference import published_models
    from grace.train import train_models

    p = prepare("athaliana", "/path/to/datasets")
    net = ensemble_network(p, published_models("athaliana"))      # the published network
    models = train_models(p.training, n_models=100)               # retrain
"""

__version__ = "1.0.0"
