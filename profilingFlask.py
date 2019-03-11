#run a profiler as a flask service

#send back info on returned hogs, linkage, tax spread

#

import profiler
from flask import jsonify

from flask import Flask, make_response
app = Flask(__name__)

@app.route('/csv/')
def download_csv():
    #return csv data of the output matrix
    csv = 'foo,bar,baz\nhai,bai,crai\n'
    response = make_response(csv)
    cd = 'attachment; filename=mycsv.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'

    return response
#build cluster heatmap json object for inchlib

#return go enrichment of top k profiles

#
