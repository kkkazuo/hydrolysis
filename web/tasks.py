import os
import sqlite3
from celery import Celery
from flask import Flask
from sklearn.externals import joblib
from mordred import Calculator, descriptors
from rdkit import Chem
import numpy as np
import pandas as pd
DATADIR = os.path.dirname(__file__)


def make_celery(app):
    celery = Celery(app.import_name, backend=app.config['CELERY_RESULT_BACKEND'],
                    broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task

    class ContextTask(TaskBase):
        abstract = True

        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)

    celery.Task = ContextTask
    return celery


flask_app = Flask(__name__)
flask_app.config.update(
    CELERY_BROKER_URL='amqp://guest@localhost//',
    CELERY_RESULT_BACKEND='amqp://guest@localhost//'
)
celery = make_celery(flask_app)


@celery.task()
def calculate(ID):
    with sqlite3.connect('main.sqlite') as conn:
        curs = conn.cursor()
        # get smiles from DB
        curs.execute("SELECT chem_index, smiles FROM result WHERE id = ?", (ID,))
        smiles_list = curs.fetchall()
        # calculate descriptors and result
        for (i, smiles) in smiles_list:
            try:
                descriptors = calculate_descriptor(smiles)
                hydrolyzability = calculate_result(descriptors)
                # update result in hydrolyzability of main.sqlite
                curs.execute(
                    "UPDATE result SET hydrolyzability = ? WHERE id = ? AND chem_index = ?",
                    (int(hydrolyzability), ID, i),
                )
            except:
                curs.execute(
                    "UPDATE result SET hydrolyzability = ? WHERE id = ? AND chem_index = ?",
                    (3, ID, i),
                )

        # update status from 0
        curs.execute("UPDATE status SET status = 1 WHERE id = ?", (ID,))


calc = Calculator(descriptors, ignore_3D=True)


# calculate descriptors with mordred
def calculate_descriptor(smiles):
    mol = Chem.MolFromSmiles(smiles)
    calculated = pd.Series(calc(mol))
    processed = calculated.apply(lambda x: np.nan if type(x) == str else float(x))
    return processed


# load imputer parameters
imp = joblib.load(os.path.join(DATADIR, './pkl_files/imp.pkl'))
# load standard scaler parameters
std = joblib.load(os.path.join(DATADIR, './pkl_files/std.pkl'))
# load model
clf = joblib.load(os.path.join(DATADIR, './pkl_files/clf.pkl'))
# load RFE
rfe = joblib.load(os.path.join(DATADIR, './pkl_files/rfe.pkl'))


# calculate result (0 or 1) with built model
def calculate_result(descriptors):
    processed = std.transform(imp.transform(descriptors))
    transformed = rfe.transform(processed)
    result = clf.predict(transformed)
    return result
