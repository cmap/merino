from flask import Flask, request, render_template
import assemble_automationify


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
            if request.form["cfg"] != '':
                raw_args.extend(["-config_filepath", request.form["cfg"]])
            if request.form["outfile"] != '':
                raw_args.extend(["-out", request.form["outfile"]])
            if request.form["pod"] != '':
                raw_args.extend(["-pod", request.form["pod"]])
            args = assemble_automationify.build_parser().parse_args(raw_args)

            try:
                assemble_automationify.main(args)
            except Exception as err:
                return "Exception: {}".format(err)
            except AssertionError as err:
                return "AssertionError: {}".format(err)

        return "Assembly Successful: GCT location {}".format(request.form["outfile"])

    return "Flask server does not handle method: {}".format(request.method)



if __name__ == '__main__':
    app.run(debug=True)