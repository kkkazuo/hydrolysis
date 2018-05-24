import csv
import io
import sqlite3
from contextlib import closing
import pandas as pd
from flask import (Flask, abort, jsonify, make_response, redirect,
                   render_template, request)
from asset.issue_id import issue_id
from tasks import calculate

app = Flask(__name__)


@app.route('/')
def home():
    return render_template('top.html')


@app.route("/file", methods=["POST"])
def post_file():
    if request.method == "POST":
        ID = issue_id()
        # get POSTED file
        # save ID and status into 'status.sqlite'
        with closing(request.files['uploaded']) as f,   sqlite3.connect('main.sqlite') as conn:
            curs = conn.cursor()
            curs.execute("INSERT INTO status (id, status) VALUES (?, ?)", (ID, 0))
            # transform file to csv format
            try:
                stream = io.StringIO(f.stream.read().decode("UTF8"), newline=None)
                csv_input = csv.reader(stream)
                # prepare for inserting data
                for ind, smi in enumerate(csv_input):
                    curs.execute(
                        'INSERT INTO result (id, chem_index, smiles) VALUES (?, ?, ?)',
                        (ID, ind, smi[0]),
                        )
                conn.commit()
                calculate.delay(ID)
                return redirect("/waiting/{}".format(ID))
            except:
                return render_template('helping.html')
    else:
        return redirect("/")


@app.route('/waiting/<ID>', methods=['GET'])
def waiting(ID):
    return render_template('waiting.html', ID=ID)


# return calculating status 0: not finished, 1: finished
@app.route('/status/<ID>')
def status(ID):
    with sqlite3.connect('main.sqlite') as conn:
        curs = conn.cursor()
        curs.execute("SELECT id, status FROM status WHERE id = ?", (ID,))
        result = curs.fetchone()
        if result is None:
            return abort(404)
        ID, status = result
        return jsonify({"status": status, "id": ID})


@app.route('/model')
def model():
    return render_template('model.html')


@app.route('/help')
def help():
    return render_template('help.html')


@app.route('/helping')
def _help():
    return render_template('helping.html')


@app.route('/result/<ID>')
def result(ID):
    # read SQL and return table
    with sqlite3.connect('main.sqlite') as conn:
        curs = conn.cursor()
        curs.execute("SELECT smiles, hydrolyzability FROM result WHERE id= ?", (ID,))
        output = curs.fetchall()
    try:
        lis = 'F,T,error'.split(',')
        output_ = [(smi, lis[i]) for (smi, i) in output]
        return render_template('result.html', result=output_, ID=ID)
    except:
        return render_template('helping.html')


@app.route('/dl/<ID>')
def download_csv(ID):
    with sqlite3.connect('main.sqlite') as conn:
        curs = conn.cursor()

        curs.execute("SELECT smiles, hydrolyzability FROM result WHERE id= ?", (ID,))
        res = curs.fetchall()
        se = pd.Series(
            data=[res[i][1] for i in range(len(res))],
            index=[res[i][0] for i in range(len(res))]
            )
        se_text = se.to_string()
    response = make_response(se_text)
    response.headers['Content-Disposition'] = 'attachment; filename=result.csv'
    response.mimetype = 'text/csv'
    return response


# for rendering
@app.route('/_footer')
def footer():
    return render_template('_footer.html')


@app.route('/_header')
def header():
    return render_template('_header.html')


if __name__ == '__main__':
    app.run(port=9999, debug=False)
