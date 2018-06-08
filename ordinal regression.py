
import mord
import  matplotlib.pyplot as plt
from sklearn import linear_model, metrics, preprocessing

from os.path import dirname, join

import numpy as np
from sklearn.datasets.base import Bunch


######

def load_housing():
    from pandas import read_csv
    """Load and return the Copenhagen housing survey dataset
       (ordinal classification).
    ==============     ==============
    Samples total                1681
    Dimensionality                  3
    Features              categorical
    Targets       ordered categorical
    ==============     ==============
    Returns
    -------
    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification targets,
        and 'DESCR', the full description of the dataset.
    Examples
    --------
    >>> from sklearn.datasets import load_housing
    >>> housing = load_housing()
    >>> print(housing.data.shape)
    (506, 13)
    """
    module_path = dirname(__file__)
    print(module_path)

    fdescr_name = 'copenhagen_housing_survey.rst'
    with open(fdescr_name) as f:
        descr_text = f.read()

    data_file_name = 'copenhagen_housing_survey.csv'
    data = read_csv(data_file_name)

    '''
    Original data set is formatted as a frequency table,
    but it's more convenient to work with the data
    as having one row per observation, below duplicates
    each obs by index based on the number the frequency ('Freq')
    of appearance
    '''
    index = np.asarray(range(0, data.shape[0])).\
        repeat(data.ix[:,'Freq'].values)
    data = data.ix[index,:]
    features = ('Infl', 'Type', 'Cont')

    return Bunch(data=data.loc[:,features],
                 target=data.loc[:,'Sat'],
                 feature_names=features,
                 DESCR=descr_text)

###
data = load_housing()
features = data.data

le = preprocessing.LabelEncoder()
le.fit(data.target)
data.target = le.transform(data.target)

features.loc[features.Infl == 'Low', 'Infl'] = 1
features.loc[features.Infl == 'Medium', 'Infl'] = 2
features.loc[features.Infl == 'High', 'Infl'] = 3

features.loc[features.Cont == 'Low', 'Cont'] = 1
features.loc[features.Cont == 'Medium', 'Cont'] = 2
features.loc[features.Cont == 'High', 'Cont'] = 3

le = preprocessing.LabelEncoder()
le.fit(features.loc[:,'Type'])
features.loc[:,'type_encoded'] = le.transform(features.loc[:,'Type'])

X, y = features.loc[:,('Infl', 'Cont', 'type_encoded')], data.target

X = X.drop(X.columns[1], axis=1)
X = X.drop(X.columns[1], axis=1)

clf1 = linear_model.LogisticRegression(
    solver='lbfgs',
    multi_class='multinomial')
clf1.fit(X, y)

print('Mean Absolute Error of LogisticRegression: %s' %
      metrics.mean_absolute_error(clf1.predict(X), y))

clf2 = mord.LogisticAT(alpha=1.)
clf2.fit(X, y)
print('Mean Absolute Error of LogisticAT %s' %
      metrics.mean_absolute_error(clf2.predict(X), y))


clf3 = mord.LogisticIT(alpha=1.)
clf3.fit(X, y)
print('Mean Absolute Error of LogisticIT %s' %
      metrics.mean_absolute_error(clf3.predict(X), y))

clf4 = mord.OrdinalRidge(alpha=1.)
clf4.fit(X, y)
print('Mean Absolute Error of LogisticSE %s' %
      metrics.mean_absolute_error(clf4.predict(X), y))
y_pred = clf4.predict(X)
plt.scatter(X, y,  color='black')
plt.plot(X, y_pred, color='blue', linewidth=3)
plt.show()