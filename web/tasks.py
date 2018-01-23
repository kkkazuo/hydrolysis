import os
import sqlite3
import json
import pandas as pd
from celery import Celery
from flask import Flask
from mordred import Calculator, descriptors
from rdkit import Chem
from rdkit.Chem import SaltRemover
from sklearn.externals import joblib

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
        for i, chemical in smiles_list:
            try:
                descriptor_values = calculate_descriptor(chemical[0])
                hydrolyzability = calculate_result(descriptor_values)

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


calc = Calculator.from_json((json.load(open("./pkl_files/calc.json", "rt"))))


# calculate descriptors with mordred
def calculate_descriptor(smiles):
    mol = Chem.MolFromSmiles(smiles)
    # remove = SaltRemover.SaltRemover()
    # removed = remove.StripMol(mol,dontRemoveEverything=True)
    calcurated = calc(mol)
    return calcurated


# load imputer parameters
imp = joblib.load(os.path.join(DATADIR, './pkl_files/Imputer.pkl'))

# load standard scaler parameters
std = joblib.load(os.path.join(DATADIR, './pkl_files/Std.pkl'))

# load model
clf = joblib.load(os.path.join(DATADIR, './pkl_files/model.pkl'))

# load RFE
rfe = joblib.load(os.path.join(DATADIR, './pkl_files/RFE.pkl'))


# calculate result (0 or 1) with built model
def calculate_result(descriptor):
    descriptor_series = pd.Series(descriptor)
    descriptor_nan = descriptor_series.apply(lambda x: '' if type(x) == str else float(x))
    processed = std.transform(imp.transform(descriptor_nan))
    transformed = rfe.transform(processed)
    result = clf.predict(transformed)
    return result
