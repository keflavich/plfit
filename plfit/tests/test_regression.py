import numpy as np

from .. import plfit

def test_issue6():
    # Issue 6 turned out to be a bug in the python 2 discrete fitter related to
    # types
    np.random.seed(0)

    test_results = [plfit(np.random.zipf(1.5,1000), discrete=True).plfit() for ii in range(20)]

    correct_results = [(11, 1.446205203835651),
                       (11, 1.5281527659073806),
                       (24, 1.4781899809035757),
                       (5, 1.4762448609017413),
                       (6, 1.5172567915706954),
                       (5, 1.4926118375988198),
                       (7, 1.4633992733678034),
                       (26, 1.4565068571094049),
                       (11, 1.463556822667798),
                       (5, 1.5208086553278282),
                       (20, 1.4684626049662981),
                       (9, 1.5119585501695973),
                       (9, 1.503821313598722),
                       (7, 1.5135408351009809),
                       (27, 1.4630831170687351),
                       (14, 1.5560600392125927),
                       (15, 1.4989900244136263),
                       (11, 1.4632893374088225),
                       (5, 1.5203141202088135),
                       (5, 1.5052229915552222)]
    
    for (xmin,alpha),(xmin_,alpha_) in zip(correct_results,test_results):
        assert xmin == xmin_
        np.testing.assert_almost_equal(alpha,alpha_)
