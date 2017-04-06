from flask import Flask, request, render_template
import assemble_automationify
import assemble
import prism_metadata
import normalization.normalize_with_prism_invariant as norm


template_directory = "/Users/elemire/Workspace/prism_pipeline"
app = Flask(__name__, template_folder=template_directory)

@app.route('/')
def make_home_page():
    return render_template("assembly_line.html")
@app.route('/', methods=["POST"])
def assembly_line():
    if request.method == "POST":
        if request.form["submit"] == "SUBMIT":
            raw_args = ["-det", request.form["det"]]
            #TODO Might be a better way to do this. As far as I could see the form always returns an empty string rather than None. So I only want to add it to the args if there was something in it.
            # Change this to a loop
            if request.form["cfg"] != '':
                raw_args.extend(["-config_filepath", request.form["cfg"]])
            if request.form["outfile"] != '':
                raw_args.extend(["-out", request.form["outfile"]])
            if request.form["pod"] != '':
                raw_args.extend(["-pod", request.form["pod"]])
            args = assemble_automationify.build_parser().parse_args(raw_args)

            try:
                outfile = assemble_automationify.main(args)
            except Exception as err:
                return "Exception: {}".format(err)
            except AssertionError as err:
                return "AssertionError: {}".format(err)

            #try:
                #norm.normalize(outfile)
            #except Exception as err:
                #return "Exception: {}".format(err)
            #except AssertionError as err:
                #return "AssertionError: {}".format(err)

        return "Assembly Successful \n GCT location {} \n NORM location {}".format(outfile, outfile[:-10] + 'NORM.gct')

    return "Flask server does not handle method: {}".format(request.method)

@app.route('/null')
def make_null_page():
    return render_template("assemble_null.html")

@app.route('/null', methods=["POST"])
def assemble_null():
    if request.method == "POST":
        if request.form["submit"] == "SUBMIT":
            assay_plates = [prism_metadata.AssayPlate(assay_plate_barcode='-666', det_plate='PNUL001_DPNULL_X1', pool_id='-P666', ignore=False)]
            assay_plates[0].det_plate_scan_time = '-666'
            null_davepool = "/Users/elemire/Workspace/flask_test/DPNULL.txt"
            null_cellset = '/Users/elemire/Workspace/flask_test/CS-1.txt'
            raw_args = ["-pmp", request.form["-pmp"],  "-prn", request.form["-prn"],
                        "-dmf", null_davepool, "-csdf", null_cellset, "-dp_csv", "DPNULL", request.form["-csv"]]
            args = assemble.build_parser().parse_args(raw_args)
            assemble.main(args, assay_plates=assay_plates)


if __name__ == '__main__':
    app.run(debug=True)
