import numpy as np

def TS_1(x1=20, samples=300, with_trend=True):
    """
    Returns:
    ========
    x,y - timestamps and values of time series for SAX exercise.
    """
    x = np.linspace(0, x1, samples)

    trend = lambda x : 0.05*x+1
    np.random.seed(5000)
    y = 0.5*np.sin(x) + np.random.uniform(size=len(x)) + (trend(x) if with_trend else 0)
    return x,y

def TS_2(x1=20, samples=300, with_trend=True):
    """
    Returns:
    ========
    x,y - timestamps and values of time series for SAX exercise.
    """
    x = np.linspace(0, x1, samples)

    trend = lambda x : 0.02*x+1
    np.random.seed(5000)
    y = 0.5*np.cos(x) + np.random.uniform(size=len(x)) + (trend(x) if with_trend else 0)
    return x,y